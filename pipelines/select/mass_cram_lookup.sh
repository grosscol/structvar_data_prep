#!/usr/bin/env bash
#  Look up a list of cram IDS
#  Extracted path generation logic from /usr/cluster/topmed/lib/perl5/Topmed/Path.pm 
# Expected Use
#   mass_cram_lookup.sh NWD100011 NWD100014 NWD100018 NWD100027 > id_paths.tsv

# Test IDS
# NWD100011 NWD100014 NWD100018 NWD100027 NWD100028 NWD100034 NWD100038 NWD100044 NWD100045 NWD100046


#############
# Functions #
#############

# Takes comma delimited list of IDS as a single string
#   "'NWD100011', 'NWD100014', 'NWD100018'"
generate_query(){
read -r -d '' TEMPLATE << HEREDOC
SELECT 
  bamid,    bamfiles.runid,    piname,        expt_sampleid,      bamname, 
  cramname, bamfiles.datayear, runs.centerid, centers.centername
FROM bamfiles
JOIN runs ON runs.runid = bamfiles.runid 
JOIN centers ON centers.centerid = runs.centerid
WHERE expt_sampleid IN 
  ($1);
HEREDOC

echo "${TEMPLATE}"
}

# Lookups for year 1 and 2 path schemes
Y1_2_TOP_DIR=('topmed10' 'topmed9' 'topmed6' 'topmed7' 'topmed9' 'topmed10/')
Y1_2_SUB_DIR=('working' 'working' 'incoming' 'incoming' 'incoming' 'incoming/')

# Arguments: bam_id, centername, piname, nwdid
# Emits path to cram file per first naming scheme for organizing files.
path_scheme_one(){
  B=$1
  ID_MOD=$((B % 2))
  TOP_DIR=${Y1_2_TOP_DIR[$ID_MOD]}
  SUB_DIR=${Y1_2_SUB_DIR[$ID_MOD]}

  echo "/net/$TOP_DIR/$SUB_DIR/mapping/results/$2/$3/b38/$4/$4.recab.cram"
}

# Arguments: bam_id, centername, piname, nwdid
# Emits path to cram file per second naming scheme for organizing files.
path_scheme_two(){
  B=$1
  ID_MOD=$((B % 6))
  TOP_DIR=${Y1_2_TOP_DIR[$ID_MOD]}
  SUB_DIR=${Y1_2_SUB_DIR[$ID_MOD]}
  echo "/net/$TOP_DIR/$SUB_DIR/mapping/results/$2/$3/b38/$4/$4.recab.cram"
}


Y3_TOP_DIR=('topmed' 'topmed2' 'topmed3' 'topmed7' 'topmed5' 'topmed6' 'topmed7' 'topmed9/')
Y3_SUB_DIR=('working' 'working' 'working' 'working' 'incoming' 'incoming' 'incoming' 'incoming/')

# Arguments: bam_id, centername, piname, nwdid
# Emits path to cram file per third naming scheme for organizing files.
path_scheme_three(){
  B=$1
  ID_MOD=$((B % 8))
  TOP_DIR=${Y3_TOP_DIR[$ID_MOD]}
  SUB_DIR=${Y3_SUB_DIR[$ID_MOD]}
  echo "/net/$TOP_DIR/$SUB_DIR/mapping/results/$2/$3/b38/$4/$4.recab.cram"
}

#############
# Main Loop #
#############

# Convert arguments (IDs) to quoted comma separated string
IDS=$(echo "'$*'" | sed "s/ /', '/g")

# Create credentials file if doesn't exist
DB_SRC_CREDS=/net/wonderland/home/grosscol/scratch/.mylogin.cnf
DB_CREDS=~/.tpmdcrdntls.cnf

cp --no-clobber ${DB_SRC_CREDS} ${DB_CREDS}
chmod 0600 ${DB_CREDS}

export MYSQL_TEST_LOGIN_FILE="${DB_CREDS}"
MYSQL_CMD_OPTS="--login-path=topmed --database=nhlbi --batch --skip-column-names"

generate_query "${IDS}" |\
  mysql ${MYSQL_CMD_OPTS} |\
  while IFS=$'\t' read -a LINEARR; do
    echo "${LINEARR[@]}"
    BAM_ID=${LINEARR[0]}
    PI_NAME=${LINEARR[2]}
    NWD_ID=${LINEARR[3]}
    CENTER_NAME=${LINEARR[8]}
    PATHS=( $(path_scheme_one $BAM_ID $CENTER_NAME $PI_NAME $NWD_ID) \
            $(path_scheme_two $BAM_ID $CENTER_NAME $PI_NAME $NWD_ID) \
            $(path_scheme_three $BAM_ID $CENTER_NAME $PI_NAME $NWD_ID) )

    echo -en "$NWD_ID\t"
    FOUND=''
    for P in "${PATHS[@]}"; do
      if stat $P 1>/dev/null 2>&1; then
        FOUND=$P
      fi
    done
    echo $FOUND
  done
