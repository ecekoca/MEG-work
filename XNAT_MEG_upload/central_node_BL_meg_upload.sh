#!/bin/bash

# Script to upload data from /megdata to xnat server using curl and XNAT REST API
# v1.1 Russell Thompson 2019-04-09
# modified by Ece K for MEG 2019


function log_msg {
# log_msg <log file> <message>

	printf "%-25s%s\n" "$(date +"%Y-%m-%d %H:%M:%S")" "${2}" >> "${1}"

}

# Define variables

# root directory - where the script lives
rdir=$(dirname $0)

# other files in root directory
subjfile="${rdir}/subject_info"
sessfile="${rdir}/BL_session_names"
logfile="${rdir}/xnat_upload_log"

# root of raw mri data path
cbupath="/megdata/cbu/ntad"


# xnat01 user credentials
# password should be kept in a file with permissions set so that only the owner
# can read it  (i.e. chmod 600 <password file>)
usr="ecekoca"
pwd=$(cat "${rdir}/xnatpwd")

# Base URL for this project on xnat01
xnat_base_url="https://central-node.dpuk.org/data"
xnat_project_url="${xnat_base_url}/archive/projects/ntad"

# Default arguments to pass to the curl command
CURL_ARGS="-Ik --silent --cookie JSESSIONID="


# check if logfile exists, if not create one. EK
if [ ! -f "${logfile}" ];then
	touch "${logfile}"
fi

# check if data exists
if [ ! -f "${subjfile}" ];then
	echo "Subject_info file not found. Quitting"
	exit
fi

# check if sessions file exists
if [ ! -f "${sessfile}" ];then
        echo "Session_names file not found. Quitting"
        exit
fi


# loop through CBU id numbers in subject info file
while IFS=$'\n' read cbuid;do


	# generate a session cookie for each new subject rather than passing username and password each time
	# this avoids the problem of having too many open sessions when uploading 100's of files

	COOKIE=$(curl -k --silent -u ${usr}:${pwd} -X POST ${xnat_base_url}/JSESSION)


	# get DPUK details for this subject as they will appear in xnat (field 2 = experiment name / xnat session name, 
	# field 4 = DPUK id number)
	# retrieve session name and NTAD ID from the subject file
	xnat_session_name=$(grep ${cbuid} "${subjfile}" | awk -F',' '{print $2}')
	xnat_subj_id=$(grep ${cbuid} "${subjfile}" | awk -F',' '{print $1}')
	cbu_subfolder=$(grep ${cbuid} "${subjfile}" | awk -F',' '{print $6}') #EK

	# get CBU meg data path for this subject
	cbusubpath=$(find "${cbupath}" -maxdepth 2 -mindepth 1 -type d -iname "${cbuid}*" 2>/dev/null)


	# if path doesn't exist, skip this subject
	if [ ! -d "${cbusubpath}" ];then
		log_msg "${logfile}" "Subject directory ${cbupath}/${cbuid}* not found. Skipping this subject."
		log_msg "${logfile}" "------------"
		continue

	# also skip if there are multiple session directories for this subject 
	elif [ $(find ${cbusubpath} -maxdepth 1 -mindepth 1 -type d -iname 20* | wc -l) -gt 1 ];then
		log_msg "${logfile}" "Subject directory ${cbusubpath} contains multiple session directories. Skipping this subject."
		log_msg "${logfile}" "------------"
         	continue

	# subject directory exists and only contains a single session
	else
		
		# URL for this subject in xnat01
		xnat_sub_url="${xnat_project_url}/subjects/${xnat_subj_id}/experiments/${xnat_session_name}"

		# test whether subject directory exist in xnat01
		curl_res=$(curl ${CURL_ARGS}${COOKIE} ${xnat_sub_url})

		# if not skip this subject
		if [ -z "$(echo "${curl_res}" | grep "HTTP/1.1 200 OK")" ];then
			log_msg "${logfile}" "Directory ${xnat_sub_url} not found. Skipping this subject"
			log_msg "${logfile}" "------------"
			continue
		else
			log_msg "${logfile}" "Processing subject ${cbuid} - DPUK id ${xnat_subj_id}"
		fi

	fi


	# loop through scans listed in session info file
	while IFS=$'\n' read sess_info;do

		xnat_scan_n=$(echo "${sess_info}" | cut -f1 -d',') # scan number used in xnat
		xnat_scan_desc=$(echo "${sess_info}" | cut -f2 -d',') # scan type / description used in xnat
		cbu_scan_desc=$(echo "${sess_info}" | cut -f3 -d',') # scan name as appears in cbu MEG data path
		xnat_url="${xnat_sub_url}/scans/${xnat_scan_n}" # xnat url for this scan

		log_msg "${logfile}" "$(printf "Scan details: %-5s%-34s%-40s" "${xnat_scan_n}" "${xnat_scan_desc}" "${cbu_scan_desc}")"
		
		# check if we have multiple scans of this type (xnat_scan_desc ends in a number)
		instance_num=$(echo "${xnat_scan_desc}" | rev | cut -f1 -d"_" | rev);
		
		# we are only expecting a single instance for this scan type
		cbuscanpath=$(find "${cbusubpath}" -maxdepth 2 -mindepth 2 -type d -iname "${cbu_scan_desc}" 2>/dev/null) #EK
		#cbuscanpath=$(find "${cbusubpath}" -maxdepth 2 -mindepth 1 -type d -iname "1*" 2>/dev/null) #EK

		log_msg "${logfile}" "cbuscanpath is ${cbuscanpath} ; cbusubpath is ${cbusubpath} ; scandescr is ${cbu_scan_desc}"
		
		# Skip this scan if can't find matching directory in cbu meg data, or if there are mutiple matching directories
		if [ "$(printf "%s\n" "${cbuscanpath}" | wc -l)" -gt 1 ];then
                        log_msg "${logfile}" "Multiple scans matching ${cbu_scan_desc} found in ${cbusubpath}. Skipping this scan"
                        log_msg "${logfile}" "------------"
                        continue
		elif [ -z "${cbuscanpath}" ] || [ ! -d "${cbuscanpath}" ];then
			log_msg "${logfile}" "Scan matching ${cbu_scan_desc} not found in ${cbuscanpath}. Skipping this scan"
			log_msg "${logfile}" "------------"
			continue
		else
			log_msg "${logfile}" "Taking data from ${cbuscanpath}"
		fi				
		

		# If it doesn't already exist, create the scan directory in xnat:
		curl_res=$(curl ${CURL_ARGS}${COOKIE} ${xnat_url})
		if [ -z "$(echo "${curl_res}" | grep "HTTP/1.1 200 OK")" ];then
			log_msg "${logfile}" "Creating directory ${xnat_url}"
			curl_res=$(curl ${CURL_ARGS}${COOKIE} -X PUT ${xnat_url}?xsiType=xnat\:megScanData\&type=${xnat_scan_desc})
			log_msg "${logfile}" "Using command: curl ${CURL_ARGS} -X PUT ${xnat_url}?xsiType=xnat\:megScanData\&type=${xnat_scan_desc}"
			
			if [ -z "$(echo "${curl_res}" | grep "HTTP/1.1 200 OK")" ];then
				log_msg "${logfile}" "Creating directory ${xnat_url} failed. Skipping this scan"
				log_msg "${logfile}" "------------"
				continue
			fi			
		else
			log_msg "${logfile}" "Directory ${xnat_url} already exists"
		fi

		# If it doesn't already exist, create the resource directory in xnat:
		curl_res=$(curl ${CURL_ARGS}${COOKIE} ${xnat_url}/resources/FIF)
		if [ -z "$(echo "${curl_res}" | grep "HTTP/1.1 200 OK")" ];then
			log_msg "${logfile}" "Creating directory ${xnat_url}/resources/FIF"

			curl_res=$(curl ${CURL_ARGS}${COOKIE} -X PUT ${xnat_url}/resources/FIF)
			log_msg "${logfile}" "curl ${CURL_ARGS} -X PUT ${xnat_url}/resources/FIF"
			
			if [ -z "$(echo "${curl_res}" | grep "HTTP/1.1 200 OK")" ];then
				log_msg "${logfile}" "Creating directory ${xnat_url}/resources/FIF failed. Skipping this scan"
				log_msg "${logfile}" "------------"
				continue
			fi
		else
			log_msg "${logfile}" "Directory ${xnat_url}/resources/FIF already exists"
		fi


		# Upload fif files
		curl_res=$(curl ${CURL_ARGS}${COOKIE} ${xnat_url}/resources/FIF)

		if [ ! -z "$(echo "${curl_res}" | grep "HTTP/1.1 200 OK")" ];then			
			log_msg "${logfile}" "Found directory ${xnat_url}/resources/FIF. Beginning data upload"
			
			# find all fif files for this scan in cbu meg data
			dcmfiles=$(find "${cbuscanpath}" -type f -iname "${cbu_scan_desc}.fif" | sort)	#EK	
			nfiles=$(printf "%s\n" ${dcmfiles} | wc -l)
			exist_count=0
			upload_count=0

			# loop through files
			for f in ${dcmfiles};do

				fname=$(basename ${f})
				
				# check if file already exists
				curl_res=$(curl ${CURL_ARGS}${COOKIE} ${xnat_url}/resources/FIF/files/${fname})


				if [ -z "$(echo "${curl_res}" | grep "HTTP/1.1 200 OK")" ];then
				# file is not there - upload from cbu meg data
					#curl -k --cookie JSESSIONID=${COOKIE} -X POST --data-binary "@${f}" "${xnat_url}/resources/FIF/files/${fname}?format=FIF&inbody=true" 
                    curl -k --cookie JSESSIONID=${COOKIE} "${xnat_url}/resources/FIF/files/${fname}" -F "${fname}"=@"${f}"
					# Check if upload has been successful
					curl_res=$(curl ${CURL_ARGS}${COOKIE} ${xnat_url}/resources/FIF/files/${fname})
					if [ ! -z "$(echo "${curl_res}" | grep "HTTP/1.1 200 OK")" ];then
						((upload_count++))		
					fi
				else
				# file already exists
					((exist_count++))				
				fi
			done
			log_msg "${logfile}" "Uploaded ${upload_count} of ${nfiles} fif files. ${exist_count} files already present"	
		fi
		log_msg "${logfile}" "------------"
	done < <(tail -n +2 "${sessfile}" | sort -n -t ',' -k1)
done < <(tail -n +2 "${subjfile}" | awk -F',' '{print $5}')
