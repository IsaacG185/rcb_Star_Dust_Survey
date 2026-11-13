# StarNet: a small CNN that learns what a star looks like on a DASCH plate cutout.
# Trained on APASS stars (known B mag) taken from the plates themselves, never on DASCH magnitudes.
# Keep this file in the project root next to rcb_paths.py. Used by 03_A_PyTorch.ipynb.

import math
import time
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
import sep
from astropy.io import fits
from scipy.spatial import cKDTree

STAMP = 41                 # size the network sees
PAD = 4                    # stored stamps are bigger so training can shift them by up to PAD px
BIG = STAMP + 2 * PAD
HALF = STAMP // 2
ASINH_S = 3.0
EDGE = 26                  # calibrators must be this far from the cutout edge
ZPAD = 30                  # zero padding around a plate so stamps never run off the edge
K_SLOTS = 52               # stamps stored per training plate: stars, artifacts, and places with no detectable star
PIX_SCALE = 1.44

CFG = dict(
    min_cal=12, win0=1.0, win_grow=1.5, min_window_stars=15,
    p_star_min=0.2, cal_off_max=3.5, cal_sig_max=0.8, cal_near_target_px=6.0,
    limit_margin=0.5, reach_n=5, reach_reject=2.0, bright_extrap=0.3, faint_extrap=0.5,
    pos_low=6.0, pos_reject=10.0, p_low=0.5,
    conf_err=0.3, conf_reach=1.0, conf_plate_gap=1.0, conf_pos=3.0, conf_p=0.7,
)


# ---------------------------------------------------------------- small helpers
class Bar:
    # Original-style ipywidgets progress bar; falls back to plain prints outside a notebook
    def __init__(self, total, label=''):
        self.total, self.label, self.t0 = max(int(total), 1), label, time.monotonic()
        self._last = 0.0
        try:
            import ipywidgets as w
            from IPython.display import display
            self.bar = w.IntProgress(value=0, min=0, max=self.total, description=label[:14])
            self.txt = w.HTML(value='')
            display(w.VBox([self.bar, self.txt]))
        except Exception:
            self.bar = self.txt = None

    def update(self, i, text=''):
        now = time.monotonic()
        if now - self._last < 0.2 and i < self.total:
            return
        self._last = now
        el = now - self.t0
        eta = el / max(i, 1) * (self.total - i)
        msg = f'{100 * i / self.total:5.1f}%  {i:,}/{self.total:,}  elapsed {el / 60:.1f} min  ETA {eta / 60:.1f} min  {text}'
        if self.bar is not None:
            self.bar.value = min(i, self.total)
            self.txt.value = f"<div style='font-family: monospace; font-size: 12px;'>{msg}</div>"
        elif i % max(self.total // 20, 1) == 0:
            print(msg)


def read_plate(path):
    with fits.open(path) as hd:
        data = np.array(hd[0].data, dtype=np.float32)
        header = hd[0].header
    if not np.isfinite(data).all():
        data = np.where(np.isfinite(data), data, np.nanmedian(data) if np.isfinite(data).any() else 0.0).astype(np.float32)
    return data, header


class Plate:
    # Background-subtracted plate in units of the local sky noise, padded with zeros
    def __init__(self, data):
        d = np.ascontiguousarray(data, dtype=np.float32)
        h, w = d.shape
        bkg = sep.Background(d, bw=min(128, w), bh=min(128, h), fw=3, fh=3)
        g = float(bkg.globalrms)
        rms = np.maximum(bkg.rms(), max(0.3 * g, 1e-3)).astype(np.float32)
        self.z = ((d - bkg.back()) / rms).astype(np.float32)
        self.h, self.w = h, w
        self.globalrms, self.globalback = g, float(bkg.globalback)
        self.zp = np.pad(self.z, ZPAD)

    def ctx(self, limit=np.nan):
        lim = 16.0 if not np.isfinite(limit) else float(limit)
        return np.array([math.log10(max(self.globalrms, 1e-3)) - 1.5,
                         math.log10(max(abs(self.globalback), 1.0)) - 3.0,
                         (lim - 15.0) / 2.0], dtype=np.float32)

    def cut(self, x, y, size=STAMP):
        # Integer-centred stamps; returns stamps (n,size,size) and the integer centres used
        x, y = np.atleast_1d(x), np.atleast_1d(y)
        xi = np.clip(np.rint(x).astype(int), -8, self.w + 8)
        yi = np.clip(np.rint(y).astype(int), -8, self.h + 8)
        hf = size // 2
        out = np.empty((len(xi), size, size), np.float32)
        for k, (a, b) in enumerate(zip(xi, yi)):
            out[k] = self.zp[b + ZPAD - hf:b + ZPAD + hf + 1, a + ZPAD - hf:a + ZPAD + hf + 1]
        return out, xi, yi


def isolated_flags(ra, dec, bmag, arcsec=12.0, dmag=2.5):
    # Not isolated if another catalogue star within arcsec is less than dmag fainter (same rule as SISaP)
    ra, dec, bmag = np.asarray(ra, float), np.asarray(dec, float), np.asarray(bmag, float)
    if len(ra) < 2:
        return np.ones(len(ra), bool)
    xy = np.c_[ra * np.cos(np.deg2rad(np.nanmean(dec))), dec] * 3600.0
    tree = cKDTree(xy)
    iso = np.ones(len(ra), bool)
    for i, nb in enumerate(tree.query_ball_point(xy, arcsec)):
        for j in nb:
            if j != i and bmag[j] < bmag[i] + dmag:
                iso[i] = False
                break
    return iso


def com_offset(z_big, r=5, thresh=3.0):
    # Centroid of the star near the middle of a BIG stamp, relative to the middle pixel; (nan, nan) if too faint
    c = BIG // 2
    win = z_big[c - r:c + r + 1, c - r:c + r + 1]
    if win.max() < 15.0:
        return np.nan, np.nan
    wgt = np.clip(win - thresh, 0, None)
    s = wgt.sum()
    if s <= 0:
        return np.nan, np.nan
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    return float((wgt * xx).sum() / s), float((wgt * yy).sum() / s)


# ---------------------------------------------------------------- training stamps for one plate
def training_plate(plate, cx, cy, cB, iso, exclude_xy=None, k_star=36, k_art=4, rng=None, min_stars=10, k_none=12):
    # One plate as fixed-size arrays, or None when it gives too few usable stars.
    # kind 1 = APASS star that is visible, 2 = artifact (trail, scratch, blob), 3 = no detectable star at the centre
    rng = rng or np.random.default_rng(0)
    h, w = plate.h, plate.w
    ok = (cx >= EDGE) & (cx < w - EDGE) & (cy >= EDGE) & (cy < h - EDGE) & iso & np.isfinite(cB)
    if exclude_xy is not None:
        ok &= np.hypot(cx - exclude_xy[0], cy - exclude_xy[1]) > 8.0
    idx = np.flatnonzero(ok)
    if len(idx) == 0:
        return None
    stamps, xi, yi = plate.cut(cx[idx], cy[idx], BIG)
    c = BIG // 2
    vis = np.array([s_[c - 4:c + 5, c - 4:c + 5].max() >= 4.0 for s_ in stamps])
    inv_idx, inv_st = idx[~vis], stamps[~vis]
    idx, stamps = idx[vis], stamps[vis]
    if len(idx) < min_stars:
        return None
    if len(idx) > k_star:                                                         # spread over the magnitude range
        order = np.argsort(cB[idx])
        pick = np.unique(np.clip((np.linspace(0, len(idx) - 1, k_star) + rng.uniform(-0.5, 0.5, k_star)).round().astype(int), 0, len(idx) - 1))
        sel = order[pick]
        idx, stamps = idx[sel], stamps[sel]
    n = len(idx)
    S = np.zeros((K_SLOTS, BIG, BIG), np.float16)
    B = np.full(K_SLOTS, np.nan, np.float32)
    kind = np.zeros(K_SLOTS, np.int8)
    C = np.zeros((K_SLOTS, 2), np.float32)
    CM = np.zeros(K_SLOTS, bool)
    S[:n] = np.arcsinh(stamps / ASINH_S)
    B[:n] = cB[idx]
    kind[:n] = 1
    for k in range(n):
        ox, oy = com_offset(stamps[k])
        if np.isfinite(ox):
            C[k], CM[k] = (ox, oy), True
    try:
        obj = sep.extract(plate.z, 3.0, err=1.0, minarea=5)
    except Exception:
        obj = []
    pos = n
    # Artifacts: detections with no catalogue star near them that are trails, scratches or blobs
    if len(obj):
        a_med = np.median(obj['a'][obj['peak'] > 15]) if (obj['peak'] > 15).sum() >= 5 else np.median(obj['a'])
        tree = cKDTree(np.c_[cx, cy])
        bad = ((obj['b'] / np.maximum(obj['a'], 1e-3)) < 0.4) | (obj['a'] > 2.5 * a_med)
        cand = np.flatnonzero(bad & (obj['peak'] > 8.0) & (obj['x'] > EDGE) & (obj['x'] < w - EDGE) & (obj['y'] > EDGE) & (obj['y'] < h - EDGE))
        cand = [i for i in cand if tree.query([obj['x'][i], obj['y'][i]])[0] > 10.0]
        if len(cand):
            cand = rng.permutation(cand)[:k_art]
            st2, _, _ = plate.cut(obj['x'][cand], obj['y'][cand], BIG)
            S[pos:pos + len(cand)] = np.arcsinh(st2 / ASINH_S)
            kind[pos:pos + len(cand)] = 2
            pos += len(cand)
    # No detectable star: catalogue stars too faint to see on this plate, plus empty sky away from every detection and catalogue star
    nn_ = 0
    if len(inv_idx):
        take = rng.permutation(len(inv_idx))[:k_none // 2]
        S[pos:pos + len(take)] = np.arcsinh(inv_st[take] / ASINH_S)
        kind[pos:pos + len(take)] = 3
        pos += len(take)
        nn_ += len(take)
    ref_xy = np.c_[np.r_[cx, obj['x'] if len(obj) else []], np.r_[cy, obj['y'] if len(obj) else []]]
    tree2 = cKDTree(ref_xy)
    tries = 0
    while nn_ < k_none and tries < 60:
        tries += 1
        rx, ry = rng.uniform(EDGE, w - EDGE), rng.uniform(EDGE, h - EDGE)
        if exclude_xy is not None and np.hypot(rx - exclude_xy[0], ry - exclude_xy[1]) < 25:
            continue
        if tree2.query([rx, ry])[0] < 14.0:
            continue
        st3, _, _ = plate.cut([rx], [ry], BIG)
        S[pos] = np.arcsinh(st3[0] / ASINH_S)
        kind[pos] = 3
        pos += 1
        nn_ += 1
    return {'S': S, 'B': B, 'kind': kind, 'C': C, 'CM': CM}


# ---------------------------------------------------------------- the network
class StarNet(nn.Module):
    def __init__(self):
        super().__init__()

        def block(i, o):
            return nn.Sequential(nn.Conv2d(i, o, 3, padding=1, bias=False), nn.BatchNorm2d(o), nn.ReLU(inplace=True))

        self.body = nn.Sequential(block(1, 24), block(24, 24), nn.MaxPool2d(2),        # 41 -> 20
                                  block(24, 48), block(48, 48), nn.MaxPool2d(2),       # 20 -> 10
                                  block(48, 96), block(96, 96), nn.MaxPool2d(2),       # 10 -> 5
                                  block(96, 128))
        self.fc = nn.Sequential(nn.Linear(128 * 5 * 5 + 3, 256), nn.ReLU(inplace=True), nn.Dropout(0.1),
                                nn.Linear(256, 64), nn.ReLU(inplace=True))
        self.head = nn.Linear(64, 5)

    def forward(self, x, ctx):
        f = self.body(x).flatten(1)
        o = self.head(self.fc(torch.cat([f, ctx], 1)))
        mu = 14.0 + o[:, 0]                       # plate-relative magnitude (the plate zero point is removed later)
        logsig = o[:, 1] - 1.0                    # log of its own error estimate
        dx, dy = 4.0 * torch.tanh(o[:, 2]), 4.0 * torch.tanh(o[:, 3])     # where the star sits relative to the stamp middle
        return mu, logsig, dx, dy, o[:, 4]        # last = star-versus-artifact logit


def _device():
    if torch.cuda.is_available():
        return torch.device('cuda')
    if getattr(torch.backends, 'mps', None) is not None and torch.backends.mps.is_available():
        return torch.device('mps')
    return torch.device('cpu')


def _batch(data, idx, rng, jitter=True):
    S = data['S'][idx]
    P, K = S.shape[:2]
    j = rng.integers(-PAD, PAD + 1, size=(P, K, 2)) if jitter else np.zeros((P, K, 2), int)
    base = np.arange(STAMP)
    rows = ((PAD + j[..., 1])[..., None] + base).reshape(P * K, STAMP)
    cols = ((PAD + j[..., 0])[..., None] + base).reshape(P * K, STAMP)
    crop = S.reshape(P * K, BIG, BIG)[np.arange(P * K)[:, None, None], rows[:, :, None], cols[:, None, :]]
    return crop.astype(np.float32), j


def _loss_and_metrics(net, data, idx, rng, device, jitter, train):
    crop, j = _batch(data, idx, rng, jitter)
    P, K = len(idx), data['S'].shape[1]
    kind = torch.from_numpy(data['kind'][idx].reshape(-1).astype(np.int64)).to(device)
    valid = kind > 0
    vidx = torch.nonzero(valid).squeeze(1)
    x = torch.from_numpy(crop).to(device)[vidx][:, None]
    ctx = torch.from_numpy(np.repeat(data['ctx'][idx], K, axis=0)).to(device)[vidx]
    mu, logsig, dx, dy, logit = net(x, ctx)

    def full(v):
        out = torch.full((P * K,), float('nan'), device=device)
        out[vidx] = v
        return out

    mu_f, ls_f, dx_f, dy_f, lg_f = full(mu), full(logsig), full(dx), full(dy), full(logit)
    B = torch.from_numpy(data['B'][idx].reshape(-1)).to(device)
    star = kind == 1
    art = kind >= 2
    r = torch.where(star, mu_f - B, torch.full_like(mu_f, float('nan'))).view(P, K)
    med = torch.nanmedian(r, dim=1, keepdim=True).values.detach()
    rc = (r - med).reshape(-1)
    sig = torch.exp(ls_f).clamp(0.03, 3.0)
    zz = rc[star] / sig[star]
    lmag = (0.5 * F.huber_loss(zz, torch.zeros_like(zz), delta=2.0, reduction='none') + torch.log(sig[star])).mean()
    C = torch.from_numpy(data['C'][idx].reshape(-1, 2)).to(device)
    cm = torch.from_numpy(data['CM'][idx].reshape(-1)).to(device) & star
    jj = torch.from_numpy(j.reshape(-1, 2).astype(np.float32)).to(device)
    if cm.any():
        loff = F.smooth_l1_loss(dx_f[cm], (C[:, 0] - jj[:, 0])[cm]) + F.smooth_l1_loss(dy_f[cm], (C[:, 1] - jj[:, 1])[cm])
    else:
        loff = torch.zeros((), device=device)
    sa = star | art
    wts = torch.where(kind[sa] == 2, torch.full_like(lg_f[sa], 3.0), torch.ones_like(lg_f[sa]))
    lbce = (F.binary_cross_entropy_with_logits(lg_f[sa], star[sa].float(), reduction='none') * wts).mean()
    loss = lmag + 0.5 * loff + 0.3 * lbce
    with torch.no_grad():
        res = rc[star].detach().cpu().numpy()
        met = {'sd': 1.4826 * np.median(np.abs(res - np.median(res))) if len(res) else np.nan,
               'off_mae': float((dx_f[cm] - (C[:, 0] - jj[:, 0])[cm]).abs().mean()) if cm.any() else np.nan,
               'acc': float(((lg_f[sa] > 0) == star[sa]).float().mean()) if sa.any() else np.nan,
               'res': res, 'B': B[star].detach().cpu().numpy()}
    return loss, met


# ---------------------------------------------------------------- the agent
class StarAgent:
    def __init__(self, device=None):
        self.device = device or _device()
        self.net = StarNet().to(self.device)
        self.meta = {}

    # -- training
    def fit(self, data, epochs=20, plates_per_batch=16, lr=2e-3, val_frac=0.05, seed=0, log=print):
        rng = np.random.default_rng(seed)
        torch.manual_seed(seed)
        n = len(data['S'])
        perm = rng.permutation(n)
        nv = max(int(n * val_frac), 1)
        val_idx, tr_idx = perm[:nv], perm[nv:]
        steps = max(len(tr_idx) // plates_per_batch, 1)
        opt = torch.optim.AdamW(self.net.parameters(), lr=lr, weight_decay=1e-4)
        sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=lr, total_steps=epochs * steps)
        rows = []
        bar = Bar(epochs * steps, 'Training')
        done = 0
        for ep in range(epochs):
            t0 = time.monotonic()
            self.net.train()
            order = rng.permutation(tr_idx)
            tl = []
            for s in range(steps):
                batch = order[s * plates_per_batch:(s + 1) * plates_per_batch]
                loss, met = _loss_and_metrics(self.net, data, batch, rng, self.device, True, True)
                opt.zero_grad()
                loss.backward()
                nn.utils.clip_grad_norm_(self.net.parameters(), 5.0)
                opt.step()
                sched.step()
                tl.append(float(loss))
                done += 1
                bar.update(done, f'epoch {ep + 1}/{epochs}  loss {np.mean(tl):.3f}')
            self.net.eval()
            vl, vsd, voff, vacc, vres, vB = [], [], [], [], [], []
            with torch.no_grad():
                for s in range(0, len(val_idx), plates_per_batch):
                    loss, met = _loss_and_metrics(self.net, data, val_idx[s:s + plates_per_batch], rng, self.device, False, False)
                    vl.append(float(loss))
                    vres.append(met['res'])
                    vB.append(met['B'])
                    voff.append(met['off_mae'])
                    vacc.append(met['acc'])
            vres, vB = np.concatenate(vres), np.concatenate(vB)
            sd = 1.4826 * np.median(np.abs(vres - np.median(vres)))
            row = {'epoch': ep + 1, 'train_loss': float(np.mean(tl)), 'val_loss': float(np.mean(vl)), 'val_sd_mag': float(sd),
                   'val_off_mae_px': float(np.nanmean(voff)), 'val_star_acc': float(np.nanmean(vacc)), 'secs': time.monotonic() - t0}
            rows.append(row)
            log(f"epoch {ep + 1:3d}  train {row['train_loss']:.3f}  val {row['val_loss']:.3f}  "
                f"val scatter {sd:.3f} mag  centre error {row['val_off_mae_px']:.2f} px  star/artifact acc {row['val_star_acc']:.2f}  {row['secs']:.0f}s")
        self.val_resid = (vres, vB)
        self.meta.update({'epochs': epochs, 'n_plates_train': int(len(tr_idx)), 'n_plates_val': int(nv)})
        return pd.DataFrame(rows)

    def save(self, path):
        torch.save({'state': self.net.state_dict(), 'meta': self.meta}, str(path))

    @classmethod
    def load(cls, path, device=None):
        ag = cls(device)
        ck = torch.load(str(path), map_location=ag.device, weights_only=False)
        ag.net.load_state_dict(ck['state'])
        ag.meta = ck.get('meta', {})
        ag.net.eval()
        return ag

    # -- inference
    @torch.no_grad()
    def predict(self, stamps, ctx, bs=512):
        # stamps: (n,STAMP,STAMP) in noise units. Returns mu, sigma, dx, dy, p_star
        self.net.eval()
        n = len(stamps)
        out = [np.empty(n, np.float32) for _ in range(5)]
        for s in range(0, n, bs):
            x = torch.from_numpy(np.arcsinh(stamps[s:s + bs] / ASINH_S).astype(np.float32))[:, None].to(self.device)
            c = torch.from_numpy(np.repeat(np.asarray(ctx, np.float32)[None], len(x), 0)).to(self.device)
            mu, ls, dx, dy, lg = self.net(x, c)
            for o, v in zip(out, (mu, torch.exp(ls).clamp(0.03, 3.0), dx, dy, torch.sigmoid(lg))):
                o[s:s + bs] = v.cpu().numpy()
        return out

    def locate(self, plate, ctx, x, y):
        # Two passes: find the star inside the stamp, re-cut on it, measure. Returns dict of arrays
        x, y = np.atleast_1d(np.asarray(x, float)), np.atleast_1d(np.asarray(y, float))
        st, xi, yi = plate.cut(x, y)
        _, _, dx1, dy1, _ = self.predict(st, ctx)
        x2, y2 = xi + dx1, yi + dy1
        st2, xj, yj = plate.cut(x2, y2)
        mu, sg, dx2, dy2, p = self.predict(st2, ctx)
        return {'mu': mu, 'sig': sg, 'p': p, 'x': xj + dx2, 'y': yj + dy2, 'off1': np.hypot(dx1, dy1)}

    def measure(self, plate, cal_x, cal_y, cal_B, tx, ty, limit=np.nan, cfg=None):
        # Magnitude of the target on one plate, calibrated on the APASS stars of that same plate
        c = {**CFG, **(cfg or {})}
        ctx = plate.ctx(limit)
        out = {'mag': np.nan, 'mag_err': np.nan, 'tier': 'reject', 'conf': '', 'flags': '', 'n_cal': 0, 'cal_rms': np.nan,
               'pos_off': np.nan, 'p_star': np.nan, 'reach_gap': np.nan, 'plate_gap': np.nan, 'mu_t': np.nan, 'sig_t': np.nan,
               'limit_mag': np.nan, 'slope': np.nan}
        cal_x, cal_y, cal_B = np.asarray(cal_x, float), np.asarray(cal_y, float), np.asarray(cal_B, float)
        keep = np.hypot(cal_x - tx, cal_y - ty) > c['cal_near_target_px']
        cal_x, cal_y, cal_B = cal_x[keep], cal_y[keep], cal_B[keep]
        if not (0 <= tx < plate.w and 0 <= ty < plate.h):
            out['flags'] = 'target_outside'
            return out
        if len(cal_x) < c['min_cal']:
            out['flags'] = 'few_calibrators'
            return out
        L = self.locate(plate, ctx, cal_x, cal_y)
        good = (L['p'] >= c['p_star_min']) & (L['off1'] <= c['cal_off_max']) & (L['sig'] <= c['cal_sig_max'])
        mu_c, B_c = L['mu'][good], cal_B[good]
        out['n_cal'] = int(good.sum())
        if good.sum() < c['min_cal']:
            out['flags'] = 'few_calibrators'
            return out
        T = self.locate(plate, ctx, [tx], [ty])
        mu_t, sig_t, p_t = float(T['mu'][0]), float(T['sig'][0]), float(T['p'][0])
        out.update({'pos_off': float(np.hypot(T['x'][0] - tx, T['y'][0] - ty)), 'p_star': p_t, 'mu_t': mu_t, 'sig_t': sig_t})
        if p_t < c['p_star_min']:                    # nothing star-like at the target: a limit at the faintest calibrator the network sees
            seen = B_c[L['p'][good] >= c['p_low']]
            out.update({'tier': 'limit', 'flags': 'not_detected', 'limit_mag': float(np.percentile(seen, 95)) if len(seen) >= 5 else np.nan})
            return out
        fit = _local_fit(mu_c, B_c, mu_t, c)
        if fit is None:
            out['flags'] = 'no_fit'
            return out
        a, b, s, n_used = fit
        order = np.sort(B_c)
        reach = order[min(c['reach_n'], len(order)) - 1]
        mag = float(a)
        err = float(max(s, sig_t) * math.sqrt(1.0 + 1.0 / max(n_used, 1)))
        flags = []
        tier = 'ok'
        if mu_t < mu_c.min() - c['bright_extrap']:
            flags.append('brighter_than_calibrators')
            tier, out['limit_mag'] = 'limit', float(order[0])
        elif mu_t > mu_c.max() + c['faint_extrap']:
            flags.append('fainter_than_calibrators')
            tier, out['limit_mag'] = 'limit', float(np.max(B_c))
        gap = float(reach - mag)
        out['reach_gap'] = gap
        if tier == 'ok' and gap > c['reach_reject']:
            flags.append('beyond_reach')
            tier = 'reject'
        if np.isfinite(limit):
            out['plate_gap'] = float(limit - mag)
            if tier == 'ok' and out['plate_gap'] < c['limit_margin']:
                flags.append('plate_limit')
                tier, out['limit_mag'] = 'limit', float(limit)
        if tier == 'ok' and out['pos_off'] > c['pos_reject']:
            flags.append('off_target')
            tier = 'reject'
        if tier == 'ok' and p_t < c['p_star_min']:
            flags.append('not_star_like')
            tier = 'reject'
        if tier == 'ok' and (out['pos_off'] > c['pos_low'] or p_t < c['p_low']):
            flags.append('low_pos_or_shape')
            tier = 'uncertain'
        out.update({'mag': mag if tier in ('ok', 'uncertain') else np.nan, 'mag_err': err, 'tier': tier, 'cal_rms': float(s), 'slope': float(b),
                    'flags': ';'.join(flags)})
        if tier in ('ok', 'uncertain'):
            pg = out['plate_gap'] if np.isfinite(out['plate_gap']) else 99.0
            high = (err <= c['conf_err'] and gap <= c['conf_reach'] and pg >= c['conf_plate_gap']
                    and out['pos_off'] <= c['conf_pos'] and p_t >= c['conf_p'])
            out['conf'] = 'high' if high else 'low'
        return out

    def detect(self, plate, ctx, nsigma=3.0, minarea=5):
        # Star finder for new plates: sep candidates, re-centred and judged by the network
        try:
            obj = sep.extract(plate.z, nsigma, err=1.0, minarea=minarea)
        except Exception:
            return pd.DataFrame()
        if len(obj) == 0:
            return pd.DataFrame()
        L = self.locate(plate, ctx, obj['x'], obj['y'])
        return pd.DataFrame({'x': L['x'], 'y': L['y'], 'mu': L['mu'], 'sigma': L['sig'], 'p_star': L['p'],
                             'peak_snr': obj['peak'], 'ellipticity': 1 - obj['b'] / np.maximum(obj['a'], 1e-3)})


def _local_fit(mu, B, mu_t, c):
    # Robust local line B = a + b*(mu - mu_t) through the calibrators near the target's mu; returns (a, b, scatter, n)
    win = c['win0']
    for _ in range(8):
        if (np.abs(mu - mu_t) <= win).sum() >= c['min_window_stars']:
            break
        win *= c['win_grow']
    m = np.abs(mu - mu_t) <= win
    if m.sum() < 6:
        return None
    x, y = mu[m] - mu_t, B[m]
    free_slope = len(x) >= 15 and (x.max() - x.min()) >= 1.0
    use = np.ones(len(x), bool)
    a, b, s = float(np.median(y - x)), 1.0, 0.3
    for _ in range(4):
        if use.sum() < 6:
            break
        if free_slope:
            A = np.c_[np.ones(use.sum()), x[use]]
            b = float(np.clip(np.linalg.lstsq(A, y[use], rcond=None)[0][1], 0.5, 1.5))
        a = float(np.median(y[use] - b * x[use]))
        res = y - (a + b * x)
        s = max(1.4826 * np.median(np.abs(res[use] - np.median(res[use]))), 0.03)
        use = np.abs(res) < 3.0 * s
    return a, b, s, int(use.sum())


# ---------------------------------------------------------------- comparing light curves
def night_pairs(jd, max_dt=0.5):
    jd = np.asarray(jd, float)
    o = np.argsort(jd)
    pairs, i = [], 0
    while i < len(o) - 1:
        if jd[o[i + 1]] - jd[o[i]] < max_dt:
            pairs.append((o[i], o[i + 1]))
            i += 2
        else:
            i += 1
    return pairs


def _rsd(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return 1.4826 * np.median(np.abs(v - np.median(v))) if len(v) else np.nan


def pull_scale(df):
    # How much the quoted errors must be multiplied so same-night pairs have a pull of 1 (independent of DASCH)
    d = df.dropna(subset=['jd', 'mag', 'mag_err']).reset_index(drop=True)
    pr = night_pairs(d['jd'])
    if len(pr) < 30:
        return 1.0, len(pr)
    a, b = np.array([p[0] for p in pr]), np.array([p[1] for p in pr])
    pull = (d['mag'].values[a] - d['mag'].values[b]) / np.hypot(d['mag_err'].values[a], d['mag_err'].values[b])
    return float(np.clip(_rsd(pull), 0.3, 3.0)), len(pr)


def compare_methods(tables, dasch=None, tol=0.5):
    # tables: {name: DataFrame(file, jd, mag, mag_err)}; dasch: DataFrame(file, jd, dasch_mag)
    # Repeatability on same-night plate pairs where every method has a point (same plates for all); agreement with DASCH on shared plates
    names = list(tables)
    base = None
    for nme in names:
        t = tables[nme][['file', 'jd', 'mag', 'mag_err']].dropna(subset=['mag']).drop_duplicates('file')
        t = t.rename(columns={'mag': f'm_{nme}', 'mag_err': f'e_{nme}'})
        base = t if base is None else base.merge(t.drop(columns='jd'), on='file')
    if dasch is not None:
        base = base.merge(dasch[['file', 'dasch_mag']].dropna().drop_duplicates('file'), on='file', how='left')
    rows = []
    pr = night_pairs(base['jd']) if len(base) else []
    for nme in names:
        t = tables[nme].dropna(subset=['mag'])
        row = {'method': nme, 'points': len(t)}
        if pr:
            a, b = np.array([p[0] for p in pr]), np.array([p[1] for p in pr])
            dd = base[f'm_{nme}'].values[a] - base[f'm_{nme}'].values[b]
            row['pairs'] = len(pr)
            row['repeat_sd'] = _rsd(dd) / math.sqrt(2)
            row['pull_sd'] = _rsd(dd / np.hypot(base[f'e_{nme}'].values[a], base[f'e_{nme}'].values[b]))
        else:
            row.update({'pairs': 0, 'repeat_sd': np.nan, 'pull_sd': np.nan})
        if dasch is not None:
            m = tables[nme].dropna(subset=['mag']).merge(dasch[['file', 'dasch_mag']].dropna(), on='file')
            dlt = m['mag'] - m['dasch_mag']
            row.update({'dasch_plates': len(m), 'within_0.5': 100 * float(np.mean(np.abs(dlt) < tol)) if len(m) else np.nan,
                        'offset': float(np.median(dlt)) if len(m) else np.nan})
        rows.append(row)
    if dasch is not None and pr:
        a, b = np.array([p[0] for p in pr]), np.array([p[1] for p in pr])
        dm = base['dasch_mag'].values
        ok = np.isfinite(dm[a]) & np.isfinite(dm[b])
        rows.append({'method': 'DASCH', 'points': int(dasch['dasch_mag'].notna().sum()), 'pairs': int(ok.sum()),
                     'repeat_sd': _rsd((dm[a] - dm[b])[ok]) / math.sqrt(2) if ok.sum() else np.nan})
    return pd.DataFrame(rows)


def merge_points(sis, ai, disagree_sigma=3.0):
    # sis, ai: DataFrame(file, jd, mag, mag_err, conf) of published points. Weighted mean where both exist.
    s = sis.rename(columns=lambda c: c + '_s' if c != 'file' else c)
    a = ai.rename(columns=lambda c: c + '_a' if c != 'file' else c)
    m = s.merge(a, on='file', how='outer')
    jd = m['jd_s'].fillna(m['jd_a'])
    both = m['mag_s'].notna() & m['mag_a'].notna()
    w_s, w_a = 1.0 / m['mag_err_s'] ** 2, 1.0 / m['mag_err_a'] ** 2
    comb = (m['mag_s'] * w_s + m['mag_a'] * w_a) / (w_s + w_a)
    ecomb = np.sqrt(1.0 / (w_s + w_a))
    emin = np.minimum(m['mag_err_s'], m['mag_err_a'])
    err = np.maximum(ecomb, 0.85 * emin)                         # the two share pixels and calibrators: no full 1/sqrt(2) gain
    dis = both & (np.abs(m['mag_s'] - m['mag_a']) > disagree_sigma * np.hypot(m['mag_err_s'], m['mag_err_a']))
    pick_s = m['mag_err_s'] <= m['mag_err_a']
    mag = np.where(both, np.where(dis, np.where(pick_s, m['mag_s'], m['mag_a']), comb), m['mag_s'].fillna(m['mag_a']))
    e = np.where(both, np.where(dis, emin, err), m['mag_err_s'].fillna(m['mag_err_a']))
    conf = np.where(both, np.where((m['conf_s'] == 'high') & (m['conf_a'] == 'high') & ~dis, 'high', 'low'),
                    m['conf_s'].fillna(m['conf_a']))
    src = np.where(both, 'both', np.where(m['mag_s'].notna(), 'sisap', 'ai'))
    return pd.DataFrame({'file': m['file'], 'jd': jd, 'mag': mag, 'mag_err': e, 'conf': conf, 'source': src,
                         'mag_sisap': m['mag_s'], 'mag_ai': m['mag_a'], 'disagree': dis}).dropna(subset=['mag', 'jd']).sort_values('jd').reset_index(drop=True)
