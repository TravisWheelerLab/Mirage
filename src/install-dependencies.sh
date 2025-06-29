#!/usr/bin/env sh


set -e


TOTAL_NUM_TESTS=27
GLOBAL_TEST_NUM=1
LOCAL_TEST_NUM=1


confirm_file_exists()
{
    FILENAME=$1
    if [ ! -e $FILENAME ];
    then
	echo "  Failed to locate file '$FILENAME'"
	exit 1
    fi
}


confirm_identical_files()
{
    FILE1=$1
    FILE2=$2
    TEST_TYPE=$3
    DIFF=$(diff $FILE1 $FILE2)
    DIFF_EXIT_CODE=$?
    if [ $DIFF_EXIT_CODE -ne 0 ];
    then
	echo "\n    ERROR ($TEST_TYPE): Files '$FILE1' and '$FILE2' do not match"
	exit 1
    fi
}


report_test_start()
{
    TEST_TYPE=$1
    printf '  [mirage2 dependencies test '
    if [ $GLOBAL_TEST_NUM -lt 10 ];
    then
	printf ' '
    fi
    printf "$GLOBAL_TEST_NUM/$TOTAL_NUM_TESTS]: $TEST_TYPE test $LOCAL_TEST_NUM..."
    GLOBAL_TEST_NUM=$(($GLOBAL_TEST_NUM + 1))
    LOCAL_TEST_NUM=$(($LOCAL_TEST_NUM + 1))
}


# sindex testing
check_sindex()
{
    FILE_NAME=$1

    INPUT_FILE=$INPUTS_DIR/$FILE_NAME
    OBSERVED_FILE=$INPUTS_DIR/$FILE_NAME.hsi
    EXPECTED_FILE=$EXP_OUTPUTS_DIR/$FILE_NAME.hsi

    if [ -e $OBSERVED_FILE ];
    then
	rm $OBSERVED_FILE
    fi

    report_test_start sindex
    $SINDEX $INPUT_FILE
    confirm_identical_files $OBSERVED_FILE $EXPECTED_FILE sindex
    echo ' passed'
}


# sstat
check_sstat()
{
    FILE_NAME=$1

    INPUT_FILE=$INPUTS_DIR/$FILE_NAME
    OBSERVED_FILE=$OBS_OUTPUTS_DIR/$FILE_NAME.sstat
    EXPECTED_FILE=$EXP_OUTPUTS_DIR/$FILE_NAME.sstat

    report_test_start sstat
    $SSTAT $INPUT_FILE > $OBSERVED_FILE
    confirm_identical_files $OBSERVED_FILE $EXPECTED_FILE sstat
    echo ' passed'
}


# sfetch
check_sfetch()
{
    RANGE_OPT=$1
    SEQ_NAME=$2
    IN_FILE_NAME=$3
    OUT_FILE_NAME=$4

    INPUT_FILE=$INPUTS_DIR/$IN_FILE_NAME
    OBSERVED_FILE=$OBS_OUTPUTS_DIR/$OUT_FILE_NAME
    EXPECTED_FILE=$EXP_OUTPUTS_DIR/$OUT_FILE_NAME

    report_test_start sfetch
    $SFETCH $RANGE_OPT $INPUT_FILE $SEQ_NAME > $OBSERVED_FILE
    confirm_identical_files $OBSERVED_FILE $EXPECTED_FILE sfetch
    echo ' passed'
}


check_spaln()
{
    AMINO_NAME=$1
    DNA_NAME=$2
    OUT_FILE_NAME=$3

    SPALN_COMP_SCRIPT=$4

    AMINO_INPUT=$INPUTS_DIR/$AMINO_NAME
    DNA_INPUT=$EXP_OUTPUTS_DIR/$DNA_NAME
    OBSERVED_FILE=$OBS_OUTPUTS_DIR/$OUT_FILE_NAME
    EXPECTED_FILE=$EXP_OUTPUTS_DIR/$OUT_FILE_NAME

    report_test_start spaln
    $SPALN -Q3 -O1 -S1 -ya3 -yz4 -yy4 $DNA_INPUT $AMINO_INPUT 1>$OBSERVED_FILE 2>/dev/null

    COMP_OUTPUT=$($SPALN_COMP_SCRIPT $EXPECTED_FILE $OBSERVED_FILE)
    COMP_EXIT_CODE=$?

    if [ $COMP_EXIT_CODE -ne 0 ];
    then
	echo "\n    ERROR (spaln): Observed hits (in '$OBSERVED_FILE') do not match expected hits (in '$EXPECTED_FILE')";
	exit 1
    fi

    echo ' passed'
}


check_blat()
{
    DNA_NAME=$1
    AMINO_NAME=$2
    OUT_FILE_NAME=$3

    DNA_INPUT=$INPUTS_DIR/$DNA_NAME
    AMINO_INPUT=$INPUTS_DIR/$AMINO_NAME
    OBSERVED_FILE=$OBS_OUTPUTS_DIR/$OUT_FILE_NAME
    EXPECTED_FILE=$EXP_OUTPUTS_DIR/$OUT_FILE_NAME

    report_test_start blat
    $BLAT -tileSize=7 -minIdentity=90 -maxIntron=1  -t=dnax -q=prot -out=blast8 -minScore=40 1>/dev/null $DNA_INPUT $AMINO_INPUT $OBSERVED_FILE
    confirm_identical_files $OBSERVED_FILE $EXPECTED_FILE blat
    echo ' passed'
}

HSI_BASE=hsi
BLAT_BASE=blat
SPALN_BASE=spaln

DEPS_DIR=dependencies

TEST_DIR=$DEPS_DIR/install-test
HSI_DIR=$DEPS_DIR/$HSI_BASE-1.0.0
BLAT_DIR=$DEPS_DIR/$BLAT_BASE
SPALN_DIR=$DEPS_DIR/${SPALN_BASE}2.2.2

TEST_TAR=$TEST_DIR.tgz
HSI_TAR=$HSI_DIR.tgz
BLAT_TAR=${BLAT_DIR}Src36.tgz
SPALN_TAR=$SPALN_DIR.tgz

confirm_file_exists $TEST_TAR
confirm_file_exists $HSI_TAR
confirm_file_exists $BLAT_TAR
confirm_file_exists $SPALN_TAR

BUILD_DIR=build

# Because we want to be sure that everything we're testing has been
# unpacked and setup in this specific environment
rm -rf $TEST_DIR
rm -rf $HSI_DIR
rm -rf $BLAT_DIR
rm -rf $SPALN_DIR

# Test input and output files
tar -xf $TEST_TAR -C $DEPS_DIR
INPUTS_DIR=$TEST_DIR/inputs
EXP_OUTPUTS_DIR=$TEST_DIR/outputs-expected
OBS_OUTPUTS_DIR=$TEST_DIR/outputs-observed

# We want to make sure we're looking at a clean slate of observations
rm -rf $OBS_OUTPUTS_DIR
mkdir $OBS_OUTPUTS_DIR

AA1_NAME=sp1.aa.fa
AA2_NAME=sp2.aa.fa
AA3_NAME=sp3.aa.fa
AA123_NAME=123.aa.fa

DNA1_NAME=sp1.dna.fa
DNA2_NAME=sp2.dna.fa
DNA3_NAME=sp3.dna.fa

GTF1=$INPUTS_DIR/sp1.gtf
GTF2=$INPUTS_DIR/sp2.gtf
GTF3=$INPUTS_DIR/sp3.gtf
SPECIES_GUIDE=$INPUTS_DIR/species-guide

confirm_file_exists $INPUTS_DIR/$AA1_NAME
confirm_file_exists $INPUTS_DIR/$AA2_NAME
confirm_file_exists $INPUTS_DIR/$AA3_NAME
confirm_file_exists $INPUTS_DIR/$AA123_NAME
confirm_file_exists $INPUTS_DIR/$DNA1_NAME
confirm_file_exists $INPUTS_DIR/$DNA2_NAME
confirm_file_exists $INPUTS_DIR/$DNA3_NAME
confirm_file_exists $GTF1
confirm_file_exists $GTF2
confirm_file_exists $GTF3
confirm_file_exists $SPECIES_GUIDE



# 1 HSI
if [[ -d $BUILD_DIR/$HSI_BASE ]]
then
    echo " HSI already built (remove $BUILD_DIR/$HSI_BASE to remake)"
else
    tar -xf $HSI_TAR -C $DEPS_DIR
    cd $HSI_DIR && cmake . && make && cd -

    SINDEX=$HSI_DIR/build/sindex
    SSTAT=$HSI_DIR/build/sstat
    SFETCH=$HSI_DIR/build/sfetch

    LOCAL_TEST_NUM=1
    check_sindex $AA1_NAME    # test 1
    check_sindex $AA2_NAME    # test 2
    check_sindex $AA3_NAME    # test 3
    check_sindex $AA123_NAME  # test 4
    check_sindex $DNA1_NAME   # test 5
    check_sindex $DNA2_NAME   # test 6
    check_sindex $DNA3_NAME   # test 7

    LOCAL_TEST_NUM=1
    check_sstat $AA1_NAME    # test  8
    check_sstat $AA2_NAME    # test  9
    check_sstat $AA3_NAME    # test 10
    check_sstat $AA123_NAME  # test 11
    check_sstat $DNA1_NAME   # test 12
    check_sstat $DNA2_NAME   # test 13
    check_sstat $DNA3_NAME   # test 14

    SP1_GENE2_DNA_NAME=sp1.g2.dna.fa
    SP1_GENE3_DNA_NAME=sp1.g3.dna.fa
    SP2_GENE2_DNA_NAME=sp2.g2.dna.fa
    SP2_GENE3_DNA_NAME=sp2.g3.dna.fa
    SP3_GENE23_DNA_NAME=sp3.g23.dna.fa

    LOCAL_TEST_NUM=1
    check_sfetch '-range 16000..38000' chr1 $DNA1_NAME $SP1_GENE2_DNA_NAME   # test 15
    check_sfetch '-range 1000..27000'  chr2 $DNA1_NAME $SP1_GENE3_DNA_NAME   # test 16
    check_sfetch '-range 13000..38000' chr1 $DNA2_NAME $SP2_GENE2_DNA_NAME   # test 17
    check_sfetch '-range 24000..1000'  chr2 $DNA2_NAME $SP2_GENE3_DNA_NAME   # test 18
    check_sfetch '-range 13000..64500' chr1 $DNA3_NAME $SP3_GENE23_DNA_NAME  # test 19

    mv -f $HSI_DIR $BUILD_DIR/$HSI_BASE
fi


# SPALN
if [[ -d $BUILD_DIR/$SPALN_BASE ]]
then
    echo " Spaln already built (remove $BUILD_DIR/$SPALN_BASE to remake)"
else
    tar -xf $SPALN_TAR -C $DEPS_DIR
    cd $SPALN_DIR/src && ./configure && make && cd -
    SPALN=$SPALN_DIR/src/spaln

    # For Spaln2.2.2 we need an comparison script that allows for
    # results to be produced in different orders.
    COMP_SCRIPT=$TEST_DIR/Compare-Spaln-Outputs.pl

    LOCAL_TEST_NUM=1
    echo '* The following 5 Spaln tests may take a minute...'
    check_spaln $AA1_NAME $SP1_GENE2_DNA_NAME sp1.g2.spaln.out   $COMP_SCRIPT # test 20
    check_spaln $AA1_NAME $SP1_GENE3_DNA_NAME sp1.g3.spaln.out   $COMP_SCRIPT # test 21
    check_spaln $AA2_NAME $SP2_GENE2_DNA_NAME sp2.g2.spaln.out   $COMP_SCRIPT # test 22
    check_spaln $AA2_NAME $SP2_GENE3_DNA_NAME sp2.g3.spaln.out   $COMP_SCRIPT # test 23
    check_spaln $AA3_NAME $SP3_GENE23_DNA_NAME sp3.g23.spaln.out $COMP_SCRIPT # test 24
    
    mv -f $SPALN_DIR $BUILD_DIR/$SPALN_BASE
fi


# BLAT
if [[ -d $BUILD_DIR/$BLAT_BASE ]]
then
    echo " Blat already built (remove $BUILD_DIR/$BLAT_BASE to remake)"
else
    tar -xzf $BLAT_TAR -C $DEPS_DIR
    cd $BLAT_DIR/lib && make && cd -
    cd $BLAT_DIR/jkOwnLib && make && cd -
    cd $BLAT_DIR/blat && mkdir -p ../bin && BINDIR=$PWD/../bin make && cd -
    BLAT=$BLAT_DIR/bin/blat
    
    LOCAL_TEST_NUM=1
    check_blat $DNA1_NAME $AA1_NAME sp1.blat.out  # test 25
    check_blat $DNA2_NAME $AA2_NAME sp2.blat.out  # test 26
    check_blat $DNA3_NAME $AA3_NAME sp3.blat.out  # test 27

    mv -f $BLAT_DIR $BUILD_DIR/$BLAT_BASE
fi


# Cleanup
rm -rf $TEST_DIR
echo "    All dependencies passed!"


# END OF SCRIPT
