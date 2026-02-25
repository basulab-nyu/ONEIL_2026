#pragma rtGlobals=1

Macro AUC_cursors()
pAUC_cursors()
proc pAUC_cursors(str_pol,str_frq,str_pls,str_out)
string str_pol
prompt str_pol "polarity (0 neg 1 pos [0])?"
string str_frq
prompt str_frq "frequency (hz [10])?"
string str_pls
prompt str_pls "pulse# ([20])?"
string str_out
prompt str_out "output?"
if(strlen(str_pol) == 0)
	str_pol = "0"
endif
if(strlen(str_frq) == 0)
	str_frq = "10"
endif
if(strlen(str_pls) == 0)
	str_pls = "20"
endif
if(strlen(str_out) == 0)
	str_out = "AUC"
endif
checkdisplayed auc_markers
if(v_flag == 1)
	removefromgraph auc_markers
endif
string str_list = sortlist(tracenamelist("",";",1),";",16)
print replacestring(";",str_list,"\r")
if(strlen(csrinfo(A)) == 0)
	cursor/p A $stringfromlist(0,str_list) 52343
endif
variable var_start = xcsr(A)
variable var_period = 1000/str2num(str_frq)
make/o/n=(itemsinlist(str_list),str2num(str_pls)) $str_out = nan
make/o/n=(str2num(str_pls),2) auc_markers = nan
variable var_index
variable var_index1
variable var_bs
var_index = 0
do
	wavestats/q/r=(var_start-100,var_start) $stringfromlist(var_index,str_list)
	var_bs = v_avg
	var_index1 = 0
	do
		duplicate/o/r=(var_start+var_index1*var_period,var_start+(var_index1+1)*var_period) $stringfromlist(var_index,str_list),temp
		temp = temp-var_bs
		if(str2num(str_pol) == 0)
			temp=0-temp
		endif
		$str_out[var_index][var_index1] = faverage(temp)
		setdimlabel 0,var_index,$stringfromlist(var_index,str_list),$str_out
		setdimlabel 1,var_index1,$("pulse_"+num2str(var_index1+1)),$str_out
		if(var_index == 0)
			auc_markers[var_index1][0] = var_start+var_index1*var_period
			auc_markers[var_index1][1] = var_bs
		endif
		var_index1+=1
	while(var_index1<str2num(str_pls))
	var_index+=1
while(var_index<itemsinlist(str_list))
appendtograph auc_markers[][1] vs auc_markers[][0]
modifygraph mode(auc_markers)=3,marker(auc_markers)=10,msize(auc_markers)=10,mrkThick(auc_markers)=2,rgb(auc_markers)=(0,0,0)
//killwaves/z temp
edit $str_out.ld
end