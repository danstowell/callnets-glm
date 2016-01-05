#!/bin/env python

###################################################
# Python script to download zebra finch call data #
###################################################

import os, errno, urllib2

# CONFIG: files to download, in the format "destination": "url"
datafiles = {
# zf4f
"data/zf4f/zcompiled_session2full.csv": "https://ndownloader.figshare.com/files/3569906",
"data/zf4f/zcompiled_session3full.csv": "https://ndownloader.figshare.com/files/3569903",
}

####################################################

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise


for destination, url in datafiles.items():
	if os.path.isfile(destination):
		print("File already exists - no need to download, so skipping: '%s'" % destination)
		continue

	# ensure folder exists
	dirname = os.path.dirname(destination)
	mkdir_p(dirname)
	# download
	u = urllib2.urlopen(url)
	block_sz = 8192
	with open(destination, 'wb') as f:
		while True:
			data = u.read(block_sz)
			if not data:
				break
			f.write(data)

