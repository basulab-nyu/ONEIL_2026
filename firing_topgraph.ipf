#pragma rtglobals=1
macro firing_topgraph()
pfiring_topgraph()
proc pfiring_topgraph(str_first,str_step,str_ovsht,str_speed,str_thr,str_bin,str_output)
string str_first
prompt str_first "first(pA)?"
string str_step
prompt str_step "step(pA)?"
string str_ovsht
prompt str_ovsht "overshoot(mV)?"
string str_speed
prompt str_speed "speed(mV/ms)?"
string str_thr
prompt str_thr "d2v/dt2(%)?"
string str_bin
prompt str_bin "bin(ms)?"
string str_output
prompt str_output "output?"
if(strlen(str_first) == 0)
	str_first = "25"
endif
if(strlen(str_step) == 0)
	str_step = "25"
endif
if(strlen(str_ovsht) == 0)
	str_ovsht = "-20"
endif
if(strlen(str_speed) == 0)
	str_speed = "5"
endif
if(strlen(str_thr) == 0)
	str_thr = "5"
endif
if(strlen(str_bin) == 0)
	str_bin = "50"
endif
killwaves/z sortwave,temp,temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10,temp11,temp12,temp13,temp14,temp15,temp16,temp17
pauseupdate
string str_list = listmatch(sortlist(tracenamelist("",";",1),";",16),"*spike_detect*")
str_list = removefromlist(greplist(str_list,"(#)"),str_list)
variable var_index = 0
if(itemsinlist(str_list)>0)
	var_index = 0
	do
		do
			checkdisplayed $stringfromlist(var_index,str_list)
			removefromgraph/z $stringfromlist(var_index,str_list)
		while(v_flag == 1)
		var_index+=1
	while(var_index<itemsinlist(str_list))
endif
resumeupdate
str_list = sortlist(tracenamelist("",";",1),";",16)
print replacestring(";",str_list,"\r")
string str_regex
string str_split
if(strlen(str_output) == 0)
	str_regex = "(w+[[:digit:]]{6}+c+[[:digit:]]{1,2}+_+[[:alpha:]]{2}+_+[[:digit:]]{3})"
	make/o/t/n=(itemsinlist(str_list)) temp = stringfromlist(p,str_list)
	var_index = 0
	do
		splitstring/e=str_regex temp[var_index],str_split
		temp[var_index] = str_split
		var_index+=1
	while(var_index<itemsinlist(str_list))
	if(dimsize(temp,0)>1)
		findduplicates/rt=temp1 temp
		str_output = temp1[0]
	else
		str_output = temp[0]
	endif
	killwaves/z temp,temp1
endif
if(strlen(csrinfo(A)) == 0 || strlen(csrinfo(B)) == 0)
	cursor/p A $stringfromlist(0,str_list) 1468
	cursor/p B $stringfromlist(0,str_list) 11468
endif
make/o/n=2 sortwave = nan
sortwave[0] = xcsr(A)
sortwave[1] = xcsr(B)
sort sortwave, sortwave
variable var_start = sortwave[0]
variable var_end = sortwave[1]
variable var_sr = 1000/dimdelta($stringfromlist(0,str_list),0)
variable var_max = 0
variable var_ahp
variable var_ahp1
variable var_ahp2
variable var_check
variable var_ss
variable var_index1
variable var_index2
var_index = 0
do
	duplicate/o/r=(var_start,var_end) $stringfromlist(var_index,str_list),temp
	findlevels/q/edge=1/d=temp1 temp,str2num(str_ovsht)
	var_max = max(var_max,dimsize(temp1,0))
	var_index+=1
while(var_index<itemsinlist(str_list))
make/o/n=(itemsinlist(str_list),var_max,16) $("spike_detect_"+str_output) = nan
setscale/p x,str2num(str_first),str2num(str_step),$("spike_detect_"+str_output)
var_index = 0
do
	setdimlabel 1,var_index,$("spike#"+num2str(var_index+1)),$("spike_detect_"+str_output)
	var_index+=1
while(var_index<dimsize($("spike_detect_"+str_output),1))
setdimlabel 2,0,thr_t,$("spike_detect_"+str_output)
setdimlabel 2,1,thr_v,$("spike_detect_"+str_output)
setdimlabel 2,2,pk_t,$("spike_detect_"+str_output)
setdimlabel 2,3,pk_v,$("spike_detect_"+str_output)
setdimlabel 2,4,ahp_t,$("spike_detect_"+str_output)
setdimlabel 2,5,ahp_v,$("spike_detect_"+str_output)
setdimlabel 2,6,hw_t1,$("spike_detect_"+str_output)
setdimlabel 2,7,hw_v1,$("spike_detect_"+str_output)
setdimlabel 2,8,hw_t2,$("spike_detect_"+str_output)
setdimlabel 2,9,hw_v2,$("spike_detect_"+str_output)
setdimlabel 2,10,amp,$("spike_detect_"+str_output)
setdimlabel 2,11,ahp,$("spike_detect_"+str_output)
setdimlabel 2,12,hw,$("spike_detect_"+str_output)
setdimlabel 2,13,isi,$("spike_detect_"+str_output)
setdimlabel 2,14,ifrq,$("spike_detect_"+str_output)
setdimlabel 2,15,lat,$("spike_detect_"+str_output)
make/o/n=(itemsinlist(str_list),max(pcsr(A),pcsr(B))-min(pcsr(A),pcsr(B))+1,2) $("spike_cutout_"+str_output) = nan
setscale/p x,str2num(str_first),str2num(str_step),$("spike_cutout_"+str_output)
setdimlabel 2,0,vm,$("spike_cutout_"+str_output)
setdimlabel 2,1,dv,$("spike_cutout_"+str_output)
make/o/n=(itemsinlist(str_list),(var_end-var_start)/str2num(str_bin),10) $("spike_bin_"+str_output) = nan
setscale/p x,str2num(str_first),str2num(str_step),$("spike_bin_"+str_output)
setscale/p y,0,str2num(str_bin),$("spike_bin_"+str_output)
var_index = 0
do
	setdimlabel 1,var_index,$("bin#"+num2str(var_index+1)),$("spike_bin_"+str_output)
	var_index+=1
while(var_index<dimsize($("spike_bin_"+str_output),1))
setdimlabel 2,0,thr_v,$("spike_bin_"+str_output)
setdimlabel 2,1,pk_v,$("spike_bin_"+str_output)
setdimlabel 2,2,pk_t,$("spike_bin_"+str_output)
setdimlabel 2,3,ahp_v,$("spike_bin_"+str_output)
setdimlabel 2,4,amp,$("spike_bin_"+str_output)
setdimlabel 2,5,ahp,$("spike_bin_"+str_output)
setdimlabel 2,6,hw,$("spike_bin_"+str_output)
setdimlabel 2,7,isi,$("spike_bin_"+str_output)
setdimlabel 2,8,ifrq,$("spike_bin_"+str_output)
setdimlabel 2,9,frq,$("spike_bin_"+str_output)
make/o/n=(1,8) $("spike_metrics_"+str_output) = nan
setdimlabel 1,0,curr,$("spike_metrics_"+str_output)
setdimlabel 1,1,thr,$("spike_metrics_"+str_output)
setdimlabel 1,2,amp,$("spike_metrics_"+str_output)
setdimlabel 1,3,ahp,$("spike_metrics_"+str_output)
setdimlabel 1,4,hw,$("spike_metrics_"+str_output)
setdimlabel 1,5,frq,$("spike_metrics_"+str_output)
setdimlabel 1,6,adapt,$("spike_metrics_"+str_output)
setdimlabel 1,7,dec,$("spike_metrics_"+str_output)
make/o/n=(1,5) $("spike_metricsRheo_"+str_output) = nan
setdimlabel 1,0,curr,$("spike_metricsRheo_"+str_output)
setdimlabel 1,1,thr,$("spike_metricsRheo_"+str_output)
setdimlabel 1,2,amp,$("spike_metricsRheo_"+str_output)
setdimlabel 1,3,ahp,$("spike_metricsRheo_"+str_output)
setdimlabel 1,4,hw,$("spike_metricsRheo_"+str_output)
var_index = 0
pauseupdate
do
	duplicate/o/r=(var_start,var_end) $stringfromlist(var_index,str_list),temp
	findlevels/q/edge=1/d=temp1 temp,str2num(str_ovsht)
	findlevels/q/edge=2/d=temp2 temp,str2num(str_ovsht)
	if(v_flag == 1)
		if(dimsize(temp1,0)!=dimsize(temp2,0))
			redimension/n=(min(dimsize(temp1,0),dimsize(temp2,0))) temp1
			redimension/n=(min(dimsize(temp1,0),dimsize(temp2,0))) temp2
		endif
		duplicate/o temp,temp3
		differentiate temp3 //temp3 = dv
		duplicate/o temp3,temp4
		differentiate temp4 //temp4 = d2v
		duplicate/o temp,temp5
		smooth 10,temp5 //temp5 = s
		duplicate/o temp3,temp6
		smooth 10,temp6 //temp6 = sdv
		duplicate/o temp4,temp7
		smooth 10,temp7 //temp7 = sd2v
		var_index1 = 0
		do
			//peak
			wavestats/q/r=(temp1[var_index1],temp2[var_index1]) temp
			$("spike_detect_"+str_output)[var_index][var_index1][%pk_t] = v_maxloc
			$("spike_detect_"+str_output)[var_index][var_index1][%pk_v] = v_max
			//ahp
			var_check = 0
			findlevel/q/edge=1/t=2/r=(temp2[var_index1],inf) temp6,0
			if(v_flag == 1)
				var_check = 1
			else
				var_ahp1 = v_levelx
			endif
			findlevel/q/edge=1/t=10/r=(temp2[var_index1],inf) temp6,0
			if(v_flag == 1)
				var_check = 1
			else
				var_ahp2 = v_levelx
			endif
			if(var_check == 1)
				$("spike_detect_"+str_output)[var_index][var_index1][] = nan
				break
			endif
			if(temp5(var_ahp1)<=temp5(var_ahp2))
				var_ahp = var_ahp1
			else
				var_ahp = var_ahp2
			endif
			$("spike_detect_"+str_output)[var_index][var_index1][%ahp_t] = var_ahp
			$("spike_detect_"+str_output)[var_index][var_index1][%ahp_v] = temp(var_ahp)
			//thr
			findlevel/q/edge=1/r=($("spike_detect_"+str_output)[var_index][var_index1][%pk_t],-inf) temp7,0
			wavestats/q/r=(v_levelx,$("spike_detect_"+str_output)[var_index][var_index1][%pk_t]) temp4
			findlevel/q/edge=1/r=($("spike_detect_"+str_output)[var_index][var_index1][%pk_t],-inf) temp7,v_max*str2num(str_thr)/100
			if(temp(v_levelx)>str2num(str_ovsht))
				findlevel/q/edge=1/r=($("spike_detect_"+str_output)[var_index][var_index1][%pk_t],-inf) temp,str2num(str_ovsht)
				findlevel/q/edge=1/r=(v_levelx,-inf) temp7,v_max*str2num(str_thr)/100
			endif
			wavestats/q/r=($("spike_detect_"+str_output)[var_index][var_index1][%pk_t]-20,$("spike_detect_"+str_output)[var_index][var_index1][%pk_t]) temp3
			if(v_max<str2num(str_speed))
				$("spike_detect_"+str_output)[var_index][var_index1][] = nan
				break
			endif
			$("spike_detect_"+str_output)[var_index][var_index1][%thr_t] = v_levelx
			$("spike_detect_"+str_output)[var_index][var_index1][%thr_v] = temp(v_levelx)
			//hw
			findlevel/q/edge=1/r=($("spike_detect_"+str_output)[var_index][var_index1][%thr_t],$("spike_detect_"+str_output)[var_index][var_index1][%pk_t]) temp,0.5*($("spike_detect_"+str_output)[var_index][var_index1][%pk_v]-$("spike_detect_"+str_output)[var_index][var_index1][%thr_v])+$("spike_detect_"+str_output)[var_index][var_index1][%thr_v]
			$("spike_detect_"+str_output)[var_index][var_index1][%hw_t1] = v_levelx
			$("spike_detect_"+str_output)[var_index][var_index1][%hw_v1] = temp(v_levelx)
			findlevel/q/edge=2/r=($("spike_detect_"+str_output)[var_index][var_index1][%pk_t],inf) temp,$("spike_detect_"+str_output)[var_index][var_index1][%hw_v1]
			$("spike_detect_"+str_output)[var_index][var_index1][%hw_t2] = v_levelx
			$("spike_detect_"+str_output)[var_index][var_index1][%hw_v2] = temp(v_levelx)
			$("spike_detect_"+str_output)[var_index][var_index1][%amp] = $("spike_detect_"+str_output)[var_index][var_index1][%pk_v]-$("spike_detect_"+str_output)[var_index][var_index1][%thr_v]
			$("spike_detect_"+str_output)[var_index][var_index1][%ahp] = $("spike_detect_"+str_output)[var_index][var_index1][%thr_v]-$("spike_detect_"+str_output)[var_index][var_index1][%ahp_v]
			$("spike_detect_"+str_output)[var_index][var_index1][%hw] = $("spike_detect_"+str_output)[var_index][var_index1][%hw_t2]-$("spike_detect_"+str_output)[var_index][var_index1][%hw_t1]
			if(var_index1>0)
				$("spike_detect_"+str_output)[var_index][var_index1][%isi] = $("spike_detect_"+str_output)[var_index][var_index1][%pk_t]-$("spike_detect_"+str_output)[var_index][var_index1-1][%pk_t]
				$("spike_detect_"+str_output)[var_index][var_index1][%ifrq] = 1000/$("spike_detect_"+str_output)[var_index][var_index1][%isi]
			endif
			$("spike_detect_"+str_output)[var_index][var_index1][%lat] = $("spike_detect_"+str_output)[var_index][var_index1][%pk_t]-var_start
			var_index1+=1
		while(var_index1<dimsize(temp1,0))
		//dvdt
		duplicate/o/r=($("spike_detect_"+str_output)[var_index][0][%thr_t],$("spike_detect_"+str_output)[var_index][dimsize(temp1,0)-1][%ahp_t]) $stringfromlist(var_index,str_list),temp8
		duplicate/o temp8,temp9
		differentiate temp9
		$("spike_cutout_"+str_output)[var_index][0,dimsize(temp8,0)-1][%vm] = temp8[q]
		$("spike_cutout_"+str_output)[var_index][0,dimsize(temp9,0)-1][%dv] = temp9[q]
		//bin
		duplicate/o/r=[var_index][][finddimlabel($("spike_detect_"+str_output),2,"lat")] $("spike_detect_"+str_output),temp10
		redimension/n=(dimsize(temp10,1)) temp10
		var_index1 = 0
		do
			extract/o/indx temp10,temp11,temp10>var_index1*str2num(str_bin) && temp10<=(var_index1+1)*str2num(str_bin)
			duplicate/o/r=[var_index][temp11[0],temp11[inf]][] $("spike_detect_"+str_output),temp12
			make/o/t/n=(dimsize($("spike_bin_"+str_output),2)) temp13 = getdimlabel($("spike_bin_"+str_output),2,p)
			var_index2 = 0
			do
				if(stringmatch(temp13[var_index2],"frq") == 1)
					wavestats/q/rmd=[][][finddimlabel(temp12,2,"pk_t")] temp12
					$("spike_bin_"+str_output)[var_index][var_index1][finddimlabel($("spike_bin_"+str_output),2,temp13[var_index2])] = 20*v_npnts/str2num(str_bin)
				else
					wavestats/q/rmd=[][][finddimlabel(temp12,2,temp13[var_index2])] temp12
					$("spike_bin_"+str_output)[var_index][var_index1][finddimlabel($("spike_bin_"+str_output),2,temp13[var_index2])] = v_avg
				endif
				var_index2+=1
			while(var_index2<dimsize(temp13,0))
			var_index1+=1
		while(var_index1<dimsize($("spike_bin_"+str_output),1))
	endif
	appendtograph $("spike_detect_"+str_output)[var_index][][%pk_v] vs $("spike_detect_"+str_output)[var_index][][%pk_t]
	appendtograph $("spike_detect_"+str_output)[var_index][][%ahp_v] vs $("spike_detect_"+str_output)[var_index][][%ahp_t]
	appendtograph $("spike_detect_"+str_output)[var_index][][%thr_v] vs $("spike_detect_"+str_output)[var_index][][%thr_t]
	appendtograph $("spike_detect_"+str_output)[var_index][][%hw_v1] vs $("spike_detect_"+str_output)[var_index][][%hw_t1]
	appendtograph $("spike_detect_"+str_output)[var_index][][%hw_v2] vs $("spike_detect_"+str_output)[var_index][][%hw_t2]
	if(var_index == 0)
		modifygraph mode($("spike_detect_"+str_output))=3,rgb($("spike_detect_"+str_output))=(65535,0,0),mrkthick($("spike_detect_"+str_output))=2
		modifygraph mode($(("spike_detect_"+str_output)+"#1"))=3,rgb($(("spike_detect_"+str_output)+"#1"))=(2,39321,1),mrkthick($(("spike_detect_"+str_output)+"#1"))=2
		modifygraph mode($(("spike_detect_"+str_output)+"#2"))=3,rgb($(("spike_detect_"+str_output)+"#2"))=(36873,14755,58982),mrkthick($(("spike_detect_"+str_output)+"#2"))=2
		modifygraph mode($(("spike_detect_"+str_output)+"#3"))=3,rgb($(("spike_detect_"+str_output)+"#3"))=(65535,43690,0),mrkthick($(("spike_detect_"+str_output)+"#3"))=2
		modifygraph mode($(("spike_detect_"+str_output)+"#4"))=3,rgb($(("spike_detect_"+str_output)+"#4"))=(65535,43690,0),mrkthick($(("spike_detect_"+str_output)+"#4"))=2
	else
		modifygraph mode($(("spike_detect_"+str_output)+"#"+num2str(var_index*5)))=3,rgb($(("spike_detect_"+str_output)+"#"+num2str(var_index*5)))=(65535,0,0),mrkthick($(("spike_detect_"+str_output)+"#"+num2str(var_index*5)))=2
		modifygraph mode($(("spike_detect_"+str_output)+"#"+num2str(var_index*5+1)))=3,rgb($(("spike_detect_"+str_output)+"#"+num2str(var_index*5+1)))=(2,39321,1),mrkthick($(("spike_detect_"+str_output)+"#"+num2str(var_index*5+1)))=2
		modifygraph mode($(("spike_detect_"+str_output)+"#"+num2str(var_index*5+2)))=3,rgb($(("spike_detect_"+str_output)+"#"+num2str(var_index*5+2)))=(36873,14755,58982),mrkthick($(("spike_detect_"+str_output)+"#"+num2str(var_index*5+2)))=2
		modifygraph mode($(("spike_detect_"+str_output)+"#"+num2str(var_index*5+3)))=3,rgb($(("spike_detect_"+str_output)+"#"+num2str(var_index*5+3)))=(65535,43690,0),mrkthick($(("spike_detect_"+str_output)+"#"+num2str(var_index*5+3)))=2
		modifygraph mode($(("spike_detect_"+str_output)+"#"+num2str(var_index*5+4)))=3,rgb($(("spike_detect_"+str_output)+"#"+num2str(var_index*5+4)))=(65535,43690,0),mrkthick($(("spike_detect_"+str_output)+"#"+num2str(var_index*5+4)))=2
	endif
	setdimlabel 0,var_index,$stringfromlist(var_index,str_list),$("spike_detect_"+str_output)
	setdimlabel 0,var_index,$stringfromlist(var_index,str_list),$("spike_cutout_"+str_output)
	setdimlabel 0,var_index,$stringfromlist(var_index,str_list),$("spike_bin_"+str_output)
	var_index+=1
while(var_index<itemsinlist(str_list))
setaxis bottom var_start-50,var_end+50
resumeupdate
//norm
duplicate/o/r=[][finddimlabel($("spike_detect_"+str_output),2,"spike#1")][finddimlabel($("spike_detect_"+str_output),2,"pk_t")] $("spike_detect_"+str_output),temp14
redimension/n=(dimsize(temp14,0)) temp14
extract/o/indx temp14,temp15,numtype(temp14) == 0
wavestats/q temp15
duplicate/o/r=[v_min,inf][][] $("spike_detect_"+str_output),$("spike_detectRheo_"+str_output)
duplicate/o/r=[v_min,inf][][] $("spike_bin_"+str_output),$("spike_binRheo_"+str_output)
setdimlabel 0,0,$getdimlabel($("spike_detect_"+str_output),0,v_min),$("spike_metricsRheo_"+str_output)
$("spike_metricsRheo_"+str_output)[0][%curr] = indextoscale($("spike_detect_"+str_output),v_min,0)
$("spike_metricsRheo_"+str_output)[0][%thr] = $("spike_detect_"+str_output)($("spike_metricsRheo_"+str_output)[0][%curr])[0][%thr_v]
$("spike_metricsRheo_"+str_output)[0][%amp] = $("spike_detect_"+str_output)($("spike_metricsRheo_"+str_output)[0][%curr])[0][%amp]
$("spike_metricsRheo_"+str_output)[0][%ahp] = $("spike_detect_"+str_output)($("spike_metricsRheo_"+str_output)[0][%curr])[0][%ahp]
$("spike_metricsRheo_"+str_output)[0][%hw] = $("spike_detect_"+str_output)($("spike_metricsRheo_"+str_output)[0][%curr])[0][%hw]
duplicate/o/r=(1.5*$("spike_metricsRheo_"+str_output)[0][%curr])[][] $("spike_detect_"+str_output),temp16
duplicate/o/r=(1.5*$("spike_metricsRheo_"+str_output)[0][%curr])[][] $("spike_bin_"+str_output),temp17
setdimlabel 0,0,$getdimlabel(temp16,0,0),$("spike_metrics_"+str_output)
$("spike_metrics_"+str_output)[0][%curr] = indextoscale(temp16,0,0)
wavestats/q/rmd=[][][finddimlabel(temp16,2,"thr_v")] temp16
$("spike_metrics_"+str_output)[0][%thr] = v_avg
wavestats/q/rmd=[][][finddimlabel(temp16,2,"amp")] temp16
$("spike_metrics_"+str_output)[0][%amp] = v_avg
wavestats/q/rmd=[][][finddimlabel(temp16,2,"ahp")] temp16
$("spike_metrics_"+str_output)[0][%ahp] = v_avg
wavestats/q/rmd=[][][finddimlabel(temp16,2,"hw")] temp16
$("spike_metrics_"+str_output)[0][%hw] = v_avg
wavestats/q/rmd=[][][finddimlabel(temp16,2,"frq")] temp16
$("spike_metrics_"+str_output)[0][%frq] = v_avg
wavestats/q/rmd=[][0.75*dimsize(temp17,1),dimsize(temp17,1)-1][finddimlabel(temp17,2,"frq")] temp17
var_ss = v_avg
wavestats/q/rmd=[][0,0.25*dimsize(temp17,1)][finddimlabel(temp17,2,"frq")] temp17
$("spike_metrics_"+str_output)[0][%adapt] = var_ss/v_avg
wavestats/q/rmd=[][0.75*dimsize(temp17,1),dimsize(temp17,1)-1][finddimlabel(temp17,2,"amp")] temp17
var_ss = v_avg
wavestats/q/rmd=[][0,0.25*dimsize(temp17,1)][finddimlabel(temp17,2,"amp")] temp17
$("spike_metrics_"+str_output)[0][%dec] = var_ss/v_avg
//FI
duplicate/o/r=[][][finddimlabel($("spike_detect_"+str_output),2,"pk_t")] $("spike_detect_"+str_output),temp18
redimension/n=(dimsize(temp18,0),dimsize(temp18,1)) temp18
make/o/n=(dimsize(temp18,0)) $("spike_fi_"+str_output) = nan
var_index = 0
do
	duplicate/o/r=[var_index][] temp18,temp19
	matrixtranspose temp19
	redimension/n=-1 temp19
	wavestats/q temp19
	$("spike_fi_"+str_output)[var_index] = v_npnts
	var_index+=1
while(var_index<dimsize(temp18,0))
setscale/p x,str2num(str_first),str2num(str_step),$("spike_fi_"+str_output)
copydimlabels/rows=0 $("spike_detect_"+str_output),$("spike_fi_"+str_output)
killwaves/z sortwave,temp,temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10,temp11,temp12,temp13,temp14,temp15,temp16,temp17,temp18,temp19
cursor/k A
cursor/k B
end
