# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
echo $argv
gs -dBATCH -dNOPAUSE -q -SDEVICE=ps2write -sOutputFile=newfile.ps $argv
