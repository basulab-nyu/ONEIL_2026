#pragma rtGlobals=1

Macro leak_access()

pleak_access()

proc pleak_access(str_output)
string str_output
prompt str_output "output waves ?"
removefromgraph/z w_base
removefromgraph/z w_peak
removefromgraph/z w_onset
removefromgraph/z w_rise_10
removefromgraph/z w_rise_90
removefromgraph/z w_rise_20
removefromgraph/z w_rise_80
removefromgraph/z w_decay_10
removefromgraph/z w_decay_90
removefromgraph/z w_decay_20
removefromgraph/z w_decay_80
removefromgraph/z w_fwhm1
removefromgraph/z w_fwhm2
removefromgraph/z w_peak_ap
removefromgraph/z w_dvdt_base
removefromgraph/z w_dvdt_peak
removefromgraph/z w_dvdt_peak_ap
string str_list = sortlist(tracenamelist("",";",1),";",16)
print replacestring(";",str_list,"\r")
make/o/n=(itemsinlist(str_list)) $("leak_"+str_output) = nan
make/o/n=(itemsinlist(str_list)) $("access_"+str_output) = nan
if(strlen(csrinfo(A)) == 0)
	cursor A $stringfromlist(0,str_list) 284.3
endif
if(strlen(csrinfo(B)) == 0)
	cursor B $stringfromlist(0,str_list) 284.4
endif
make/o/n=2 sortwave = nan
sortwave[0] = xcsr(A)
sortwave[1] = xcsr(B)
sort sortwave, sortwave
variable var_index = 0
do
	$("leak_"+str_output)[var_index] = $stringfromlist(var_index,str_list)(sortwave[0])
	$("access_"+str_output)[var_index] = $stringfromlist(var_index,str_list)(sortwave[0]) - $stringfromlist(var_index,str_list)(sortwave[1])
	setdimlabel 0,var_index,$stringfromlist(var_index,str_list),$("leak_"+str_output)
	setdimlabel 0,var_index,$stringfromlist(var_index,str_list),$("access_"+str_output)
	var_index+=1
while(var_index<itemsinlist(str_list))
end