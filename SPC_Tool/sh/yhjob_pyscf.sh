

cont=0
num_freq=216
num_mission=10
name="DyCp"
shopt -s extglob

for ((imod=1;imod<=$num_freq;imod++));do
{
	cd $imod/
	tempFile="6/${name}_rohf.chk"
	if [ ! -f "$tempFile" ];then
		echo $PWD
		#cp ../eq/*_rohf.chk ./6/
		#cp ../eq/*_rohf.chk ./5/
		qsub -N $name$imod pyjob.sh  
		#sbatch -J $name$imod pyjob.sh
		cont=$[cont+1]
	else
		echo "Have submited mode $imod"
	fi
	cd ../
	if [ $cont -eq $num_mission ];then
		break
	fi
}
done
