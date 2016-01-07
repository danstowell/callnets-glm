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
# gill et al
"data/gilletal_trial2/gill_rawfile_1ind.csv":  "https://ndownloader.figshare.com/files/3662322?private_link=73cbdf96ab156a0f0a69",
"data/gilletal_trial2/gill_rawfile_1type.csv": "https://ndownloader.figshare.com/files/3662325?private_link=73cbdf96ab156a0f0a69",
"data/gilletal_trial2/gill_rawfile_4ind.csv":  "https://ndownloader.figshare.com/files/3662328?private_link=73cbdf96ab156a0f0a69",
"data/gilletal_trial2/gill_rawfile_4type.csv": "https://ndownloader.figshare.com/files/3662331?private_link=73cbdf96ab156a0f0a69",
"data/gilletal_trial2/gill_rawfile_5ind.csv":  "https://ndownloader.figshare.com/files/3662334?private_link=73cbdf96ab156a0f0a69",
"data/gilletal_trial2/gill_rawfile_5type.csv": "https://ndownloader.figshare.com/files/3662337?private_link=73cbdf96ab156a0f0a69",
"data/gilletal_trial2/gill_rawfile_6ind.csv":  "https://ndownloader.figshare.com/files/3662340?private_link=73cbdf96ab156a0f0a69",
"data/gilletal_trial2/gill_rawfile_6type.csv": "https://ndownloader.figshare.com/files/3662343?private_link=73cbdf96ab156a0f0a69",
"data/gilletal_trial2/key_gill_ind.csv":       "https://ndownloader.figshare.com/files/3662346?private_link=73cbdf96ab156a0f0a69",
"data/gilletal_trial2/key_gill_type.csv":      "https://ndownloader.figshare.com/files/3662349?private_link=73cbdf96ab156a0f0a69",
"data/gilletal_trial2/README.txt":             "https://ndownloader.figshare.com/files/3662352?private_link=73cbdf96ab156a0f0a69",
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
	print("Downloading: '%s'" % destination)
	u = urllib2.urlopen(url)
	block_sz = 8192
	with open(destination, 'wb') as f:
		while True:
			data = u.read(block_sz)
			if not data:
				break
			f.write(data)

