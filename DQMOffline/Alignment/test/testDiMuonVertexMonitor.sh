#!/bin/bash

 function die { echo $1: status $2; exit $2; }

if [ "${SCRAM_TEST_NAME}" != "" ] ; then
  mkdir ${SCRAM_TEST_NAME}
  cd ${SCRAM_TEST_NAME}
fi

echo -e " Tesing on Z->mm \n\n"

cmsRun ${LOCAL_TEST_DIR}/DiMuonTkAlDQMValidator_cfg.py resonance=Z  || die "Failure using DiMuonTkAlDQMValidator_cfg.py resonance=Z" $?
cmsRun ${LOCAL_TEST_DIR}/DiMuonTkAlDQMHarvester_cfg.py resonance=Z || die "Failure using DiMuonTkAlDQMHarvester_cfg.py resonance=Z" $? 

echo -e " Testing on J/psi -> mm \n\n"

cmsRun ${LOCAL_TEST_DIR}/DiMuonTkAlDQMValidator_cfg.py resonance=Jpsi || die "Failure using DiMuonTkAlDQMValidator_cfg.py resonance=Jpsi" $?
cmsRun ${LOCAL_TEST_DIR}/DiMuonTkAlDQMHarvester_cfg.py resonance=Jpsi || die "Failure using DiMuonTkAlDQMHarvester_cfg.py resonance=Jpsi" $? 

echo -e " Testing on Upsilon -> mm \n\n"

cmsRun ${LOCAL_TEST_DIR}/DiMuonTkAlDQMValidator_cfg.py resonance=Upsilon || die "Failure using DiMuonTkAlDQMValidator_cfg.py resonance=Upsilon" $?
cmsRun ${LOCAL_TEST_DIR}/DiMuonTkAlDQMHarvester_cfg.py resonance=Upsilon || die "Failure using DiMuonTkAlDQMHarvester_cfg.py resonance=Upsilon" $? 
