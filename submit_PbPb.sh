if [ $# -ne 0 ]
then
  echo "Usage: ./psort.sh <trackqual> <file-list> <tag> <nmin> <nmax> <pttrigmin> <pttrigmax> <ptassmin> <ptassmax>"
  exit 1
fi

now="RAA_$(date +"%Y_%m_%d__%H_%M_%S")"
njobs=700

mkdir $now
cp fileLists/PbPb_fileList.txt $now/fileList.txt 
cp fileLists/testPbPb_fileList.txt $now/testFileList.txt
cp -r TrkCorr_* $now
tar -cvzf $now/Corrections.tar.gz TrkCorr_*
cp run.sh $now

cat run.condor | sed "s@log_flag@$now@g" | sed "s@dir_flag@$PWD/$now@g" | sed "s@user_flag@$USER@g" |  sed "s@arglist@$njobs 0@g" | sed "s@transfer_filelist@run.exe@g" | sed "s@njobs@$njobs@g" > $now/run.condor

NAME="countTracks.C"
g++ $NAME $(root-config --cflags --libs)  -Wall -O2 -o "run.exe"
cp run.exe $now
rm run.exe
echo finished compilation
echo
cat $now/run.condor
echo 
echo condor_submit $now/run.condor
echo
# condor_submit $now/run.condor

