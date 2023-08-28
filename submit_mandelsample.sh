#!/bin/sh                                                                                                                                                                                             

echo "Executing MandelSample Macro"
echo "____________________________"
MyDataInDir=/raid12/zbaldwin/MandelSample

COUNTER=1
MyRun_Acc=(${MyDataInDir}/pi0eta_pi0pi0pippim_MCRecon_FlatTree_noPhotos_3*_SetThrownArray.root)
MyRun_Gen=(${MyDataInDir}/pi0eta_pi0eta_pi0pi0pippim_MCThrown_Tree_noPhotos_3*_SetThrownArray.root)

count=${#MyRun_Acc[@]}
echo $count

#FOR ALL THE RUN NUMBERS
for ((i=0; i<${#MyRun_Acc[@]} ;i++)); do

echo ${MyRun_Acc[i]}
echo ${MyRun_Gen[i]}
echo ________________________________________________________________________

python3 mandelsample.py pi0eta_Data_FlatTree.root ${MyRun_Acc[i]} ${MyRun_Gen[i]}

echo ______________________
printf "Files ran: %d\n" $COUNTER
let COUNTER++
echo ______________________

done


