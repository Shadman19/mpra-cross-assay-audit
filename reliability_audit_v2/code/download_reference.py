"""Download the public hg38 two-bit reference with the publisher checksum."""
from common import ROOT,fetch_checked
if __name__=='__main__':
    fetch_checked('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit',ROOT/'runtime/hg38.2bit','dcc3ea27079aa6dc3f9deccd7275e0f8',algorithm='md5')
    print('hg38.2bit matches the UCSC publisher MD5.')
