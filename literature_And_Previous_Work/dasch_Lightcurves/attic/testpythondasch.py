# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

import dasch
print dasch.GetMaxGSCBin.__doc__
print dasch.GetMaxGSCBin()
print dasch.GetGSCBin.__doc__
print dasch.GetGSCBin(0.0,89.99)
print dasch.GetGSCBin(359.9,89.99)
print dasch.GetBinCenter.__doc__
print dasch.GetBinCenter(52135)
print dasch.InitDatabaseAccess.__doc__
print dasch.InitDatabaseAccess("apass")
print dasch.GetFileSummaryMagnitudes.__doc__
measurements = dasch.GetFileSummaryMagnitudes(147305112)
print "GetFileSummaryMagnitudes found ",measurements," measurements"
print dasch.GetFileStarImage.__doc__
count = 0
while (count < measurements) & (count < 10):
    image = dasch.GetFileStarImage(count)
    print "image number ",count," REF is ",image.get("REF")," ra: ",image.get("ra")," dec: ",image.get("dec")
    count += 1

print image.items()
    
        

