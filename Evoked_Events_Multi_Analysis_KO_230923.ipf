#pragma rtGlobals=1
#include <Multi-peak Fitting 2.0>

Macro Evoked_Events_Multi_Analysis_KO_230923()

pEEMA_KO()

proc pEEMA_KO(str_base_w,str_peak_w,str_pol,str_thr,str_multi,str_pulse,str_frq,str_dvdt,str_output)
//declaring parameters
string str_base_w
prompt str_base_w "baseline window (ms) ?"
string str_peak_w
prompt str_peak_w "peak window (ms) ?"
string str_pol
prompt str_pol "polarity (positive [1] / negative [0]) ?"
string str_thr
prompt str_thr "spike detection (threshold [level] / override [NaN]) ?"
string str_multi
prompt str_multi "multipeak detection (smoothing [factor] / override [NaN]) ?"
string str_pulse
prompt str_pulse "multiple stims (pulses [number]) ?"
string str_frq
prompt str_frq "stim frequency (multi pulse [hz] / single pulse [NaN]) ?"
string str_dvdt
prompt str_dvdt "dV/dt analysis (threshold [level] / override [NaN]) ?"
string str_output
prompt str_output "output waves ?"
//initializing tiling
string/g str_graph_list = sortlist(winlist("*",";","WIN:1"),";",16)
string/g str_table_list = sortlist(winlist("*",";","WIN:2"),";",16)
variable/g var_graph_last = itemsinlist(str_graph_list)-1
string/g str_graph_last = winname(0,1)
variable/g var_table_last = itemsinlist(str_table_list)-1
string/g str_table_last = winname(0,2)
getwindow $str_graph_last, wsize
variable/g var_left_g = v_left + 5.25
variable/g var_top_g = v_top + 22.5
variable/g var_right_g = v_right + 5.25
variable/g var_bottom_g = v_bottom + 22.5
getwindow $str_table_last, wsize
variable/g var_left_t = v_left + 5.25
variable/g var_top_t = v_top + 22.5
variable/g var_right_t = v_right + 5.25
variable/g var_bottom_t = v_bottom + 22.5
//cleaning up
silent 1
pauseupdate
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
if(strlen(listmatch(tracenamelist("",";",1),"w_multi_pk_v_*"))>0)
	variable/g var_index = 0
	do
		removefromgraph/z $stringfromlist(var_index,wavelist("w_multi_pk_v_*",";",""))
		var_index+=1
	while(var_index<itemsinlist(wavelist("w_multi_pk_v_*",";","")))
endif
if(strlen(listmatch(tracenamelist("",";",1),"w*dvdt*"))>0)
	variable/g var_index = 0
	do
		removefromgraph/z $stringfromlist(var_index,wavelist("w*dvdt*",";",""))
		var_index+=1
	while(var_index<itemsinlist(wavelist("w*dvdt*",";","")))
endif
modifygraph axisenab(left)={0,1}
//listing input traces
string/g str_input_list = sortlist(tracenamelist("",";",1),";",16)
//preparing output waves
string/g str_output_list = "bs_"+str_output+";"+"bsT_"+str_output+";"+"amp_"+str_output+";"+"pk_"+str_output+";"+"pkT_"+str_output+";"+"APpb_"+str_output+";"+"APt_"+str_output+";"+"ostB_"+str_output+";"+"ostR_"+str_output+";"+"lagB_"+str_output+";"+"lagR_"+str_output+";"+"rs19_"+str_output+";"+"rs28_"+str_output+";"+"dc19_"+str_output+";"+"dc28_"+str_output+";"+"ar19_"+str_output+";"+"hw_"+str_output+";"+"Npk_"+str_output+";"+"sum_"+str_output
string/g str_display_list = "w_base;w_base_t;w_peak;w_peak_t;w_onset;w_onset_t;w_rise_10;w_rise_10_t;w_rise_90;w_rise_90_t;w_rise_20;w_rise_20_t;w_rise_80;w_rise_80_t;w_decay_10;w_decay_10_t;w_decay_90;w_decay_90_t;w_decay_20;w_decay_20_t;w_decay_80;w_decay_80_t;w_fwhm1;w_fwhm1_t;w_fwhm2;w_fwhm2_t;w_peak_ap;w_peak_ap_t;"
string/g str_multi_list = "multi_amp_"+str_output+";"+"multi_pk_"+str_output+";"+"multi_pkT_"+str_output
string/g str_dvdt_list = "dvdt_bs_"+str_output+";"+"dvdt_bsT_"+str_output+";"+"dvdt_amp_"+str_output+";"+"dvdt_pk_"+str_output+";"+"dvdt_pkT_"+str_output+";"+"dvdt_APpb_"+str_output+";"+"dvdt_APt_"+str_output+";"+"dvdt_APamp_"+str_output
string/g str_display_dvdt_list = "w_dvdt_base;w_dvdt_base_t;w_dvdt_peak;w_dvdt_peak_t;w_dvdt_peak_ap;w_dvdt_peak_ap_t;"
//printing I/O info
print("\rdetection parameters:")
print ("baseline averaging window = "+str_base_w+"ms")
string/g str_bs_label = str_base_w+"ms"
print ("peak averaging window = "+str_peak_w+"ms")
string/g str_pk_label = str_base_w+"ms"
if (str2num(str_pol) == 1)
	print ("polarity = positive")
	string/g str_pol_label = "positive"
else
	print ("polarity = negative")
	string/g str_pol_label = "negative"
endif
if(numtype(str2num(str_thr)) == 2)
	print ("spike threshold = overriden")
	string/g str_thr_label = "overriden"
else
	print ("spike threshold = "+str_thr)
	string/g str_thr_label = str_thr
endif
if(numtype(str2num(str_multi)) == 2)
	print ("multipeak detection factor = overriden")
	string/g str_multi_label = "overriden"
else
	print ("multipeak smoothing factor = "+str_multi)
	string/g str_multi_label = str_multi
endif
if(numtype(str2num(str_pulse)) == 2)
	str_pulse = "1"
endif
print ("number of pulses = "+str_pulse)
string/g str_pls_label = str_pulse+" pulses"
if(numtype(str2num(str_frq)) == 2)
	print ("stim frequency = N/A")
	string/g str_frq_label = "N/A"
else
	print ("stim frequency = "+str_frq+"hz")
	string/g str_frq_label = str_frq+"hz"
endif
if(numtype(str2num(str_dvdt)) == 2)
	print ("dV/dt analysis = overriden")
	string/g str_dvdt_label = "overriden"
else
	print ("dV/dt threshold = "+str_dvdt)
	string/g str_dvdt_label = str_dvdt
endif
print "\rinput traces:"
print replacestring(";",str_input_list,"\r")
print "\routput waves:"
print replacestring(";",str_output_list,"\r")
print ("data_array_"+str_output)
if(numtype(str2num(str_multi)) == 0)
	print replacestring(";",str_multi_list,"\r")
endif
if(numtype(str2num(str_dvdt)) == 0)
	print replacestring(";",str_dvdt_list,"\r")
endif
//grabbing marker position
make/o/n=4 sortwave = NaN
sortwave[0] = xcsr(A)
sortwave[1] = xcsr(B)
sortwave[2] = xcsr(C)
sortwave[3] = xcsr(D)
sort sortwave, sortwave
variable/g var_stim = sortwave[0]
variable/g var_start = sortwave[1]
variable/g var_search = sortwave[2]
variable/g var_end = sortwave[3]
killwaves/z sortwave
//setting up
string/g str_params = replacestring(",","str_bs_label,str_pk_label,str_pol_label,str_thr_label,str_multi_label,str_pls_label,str_frq_label,str_dvdt_label,str_output",";")
string/g str_labels = "baseline averaging window;peak averaging window;polarity;spike threshold;multipeak detection factor;number of pulses;stim frequency;dV/dt threshold;output name"
make/o/t/n=(max(itemsinlist(str_params),itemsinlist(str_input_list),itemsinlist(str_output_list)),3) $("params_"+str_output)
variable/g var_scroll = 0
do
	$("params_"+str_output)[var_scroll][0] = $stringfromlist(var_scroll,str_params)
	setdimlabel 0,var_scroll,$stringfromlist(var_scroll,str_labels),$("params_"+str_output)
	var_scroll += 1
while (var_scroll < itemsinlist(str_params))
variable/g var_scroll = 0
do
	$("params_"+str_output)[var_scroll][1] = stringfromlist(var_scroll,str_input_list)
	var_scroll += 1
while (var_scroll < itemsinlist(str_input_list))
variable/g var_scroll = 0
do
	$("params_"+str_output)[var_scroll][2] = stringfromlist(var_scroll,str_output_list)
	var_scroll += 1
while (var_scroll < itemsinlist(str_output_list))
setdimlabel 1,0,$"detection parameters",$("params_"+str_output)
setdimlabel 1,1,$"input traces",$("params_"+str_output)
setdimlabel 1,2,$"output waves",$("params_"+str_output)
variable/g var_base_w = str2num(str_base_w)
variable/g var_peak_w = str2num(str_peak_w)
variable/g var_pol = str2num(str_pol)
variable/g var_thr = str2num(str_thr)
variable/g var_multi = str2num(str_multi)
variable/g var_dvdt = str2num(str_dvdt)
if (numtype(str2num(str_frq)) == 2)
	variable/g var_frq = 1
else
	variable/g var_frq = 1000/str2num(str_frq)
endif
variable/g var_pulse = str2num(str_pulse)
variable/g var_total = itemsinlist(str_input_list)
if (numtype(var_thr) == 2 && var_pol == 0)
	var_thr = -inf
endif
if (numtype(var_thr) == 2 && var_pol == 1)
	var_thr = inf
endif
variable/g var_scroll = 0
do
	make/o/n=(var_total,var_pulse) $stringfromlist(var_scroll,str_output_list) = NaN
	var_scroll += 1
while (var_scroll < itemsinlist(str_output_list))
make/o/n=(var_total,var_pulse,itemsinlist(str_output_list)) $("data_array_"+str_output)
if(numtype(str2num(str_dvdt)) == 0)
	variable/g var_scroll = 0
	do
		make/o/n=(var_total,var_pulse) $stringfromlist(var_scroll,str_dvdt_list) = NaN
		var_scroll += 1
	while (var_scroll < itemsinlist(str_dvdt_list))
endif
if(numtype(str2num(str_multi)) == 0)
	variable/g var_scroll = 0
	do
		make/o/n=(var_total,var_pulse,1000) $stringfromlist(var_scroll,str_multi_list) = NaN
		var_scroll += 1
	while (var_scroll < itemsinlist(str_multi_list))
endif
variable/g var_scroll = 0
do
	make/o/n=(var_total,var_pulse) $stringfromlist(var_scroll,str_display_list) = NaN
	var_scroll += 1
while (var_scroll < itemsinlist(str_display_list))
make/o/n=1000 w_colors = x
if((numtype(str2num(str_dvdt)) == 0))
	make/o/n=(var_total,var_pulse) w_dvdt_base = nan
	make/o/n=(var_total,var_pulse) w_dvdt_base_t = nan
	make/o/n=(var_total,var_pulse) w_dvdt_peak = nan
	make/o/n=(var_total,var_pulse) w_dvdt_peak_t = nan
	make/o/n=(var_total,var_pulse) w_dvdt_peak_ap = nan
	make/o/n=(var_total,var_pulse) w_dvdt_peak_ap_t = nan
endif
//analyzing evoked events
variable/g var_index = 0
do //looping over sweeps
	variable/g var_fail = 0
	string/g str_waves = stringfromlist(var_index,str_input_list)
	variable/g var_plscount = 0
	do //looping over pulses
//baseline
		wavestats/q/r=(var_plscount*var_frq+var_start-0.5*var_base_w,var_plscount*var_frq+var_start+0.5*var_base_w) $str_waves
		variable/g var_bs = v_avg
		$("bs_"+str_output)[var_index][var_plscount] = v_avg
		$("bsT_"+str_output)[var_index][var_plscount] = var_plscount*var_frq+var_start
		w_base[var_index][var_plscount] = $str_waves(var_plscount*var_frq+var_start)
		w_base_t[var_index][var_plscount] = var_plscount*var_frq+var_start
//finding local peak
		wavestats/q/r=(var_plscount*var_frq+var_start,var_plscount*var_frq+var_search) $str_waves
		if (var_pol == 1) //upward
			variable/g var_peak_t = v_maxloc
			variable/g var_peak = v_max
		else //downward
			variable/g var_peak_t = v_minloc
			variable/g var_peak = v_min
		endif
		if (round(var_peak_t - (var_plscount*var_frq+var_start))<1) //failsafe
			var_peak_t = var_peak_t+1
			variable/g var_fail = 1
		endif
		if ((var_pol == 1 && var_peak-var_bs < var_thr) || (var_pol == 0 && var_peak-var_bs > var_thr)) //subthreshold
//amplitude
			wavestats/q/r=(var_peak_t-0.5*var_peak_w,var_peak_t+0.5*var_peak_w) $str_waves
			$("amp_"+str_output)[var_index][var_plscount] = v_avg-$("bs_"+str_output)[var_index][var_plscount]
			$("pk_"+str_output)[var_index][var_plscount] = v_avg
			$("APpb_"+str_output)[var_index][var_plscount] = 0
			$("pkT_"+str_output)[var_index][var_plscount] = var_peak_t-(var_plscount*var_frq+var_stim)
			if (var_plscount ==0) //if first peak
				$("sum_"+str_output)[var_index][var_plscount] = $("amp_"+str_output)[var_index][var_plscount] //variable is just the amplitude
			else
				$("sum_"+str_output)[var_index][var_plscount] =$("pk_"+str_output)[var_index][var_plscount] -$("bs_"+str_output)[var_index][0] //subtract baseline from pulse 1 
			endif
			
			if (var_pol == 1) //upward
				w_peak_t[var_index][var_plscount] = v_maxloc
				w_peak[var_index][var_plscount] = v_max
			else //downward
				w_peak_t[var_index][var_plscount] = v_minloc
				w_peak[var_index][var_plscount] = v_min
			endif
//kinetics
			edgestats/q/r=(var_plscount*var_frq+var_start,var_peak_t) $str_waves
			if (V_flag == 2) //failsafe
				variable/g var_fail = 1
			else //rise1090
				variable/g var_r10 = v_edgeloc1
				$("rs19_"+str_output)[var_index][var_plscount] = v_edgeloc3-v_edgeloc1
				$("ostR_"+str_output)[var_index][var_plscount] = v_edgeloc1
				$("lagR_"+str_output)[var_index][var_plscount] = v_edgeloc1-(var_plscount*var_frq+var_stim)
				w_rise_10[var_index][var_plscount] = v_edgelvl1
				w_rise_10_t[var_index][var_plscount] = v_edgeloc1
				w_rise_90[var_index][var_plscount] = v_edgelvl3
				w_rise_90_t[var_index][var_plscount] = v_edgeloc3
			endif
			edgestats/f=0.2/q/r=(var_plscount*var_frq+var_start,var_peak_t) $str_waves
			if (V_flag == 2) //failsafe
				variable/g var_fail = 1
			else //rise2080
				$("rs28_"+str_output)[var_index][var_plscount] = v_edgeloc3-v_edgeloc1
				w_rise_20[var_index][var_plscount] = v_edgelvl1
				w_rise_20_t[var_index][var_plscount] = v_edgeloc1
				w_rise_80[var_index][var_plscount] = v_edgelvl3
				w_rise_80_t[var_index][var_plscount] = v_edgeloc3
			endif
			findlevel/q/r=(var_peak_t,var_plscount*var_frq+var_start) $str_waves,$("bs_"+str_output)[var_index][var_plscount] //onset
			$("ostB_"+str_output)[var_index][var_plscount] = v_levelx
			$("lagB_"+str_output)[var_index][var_plscount] = v_levelx-(var_plscount*var_frq+var_stim)
			w_onset[var_index][var_plscount] = $str_waves(v_levelx)
			w_onset_t[var_index][var_plscount] = v_levelx
			edgestats/q/r=(var_peak_t,var_plscount*var_frq+var_end) $str_waves
			if (V_flag == 2) //failsafe
				variable/g var_fail = 1
			else //decay1090
				variable/g var_d90 = v_edgeloc3
				$("dc19_"+str_output)[var_index][var_plscount] = v_edgeloc3-v_edgeloc1
				w_decay_10[var_index][var_plscount] = v_edgelvl1
				w_decay_10_t[var_index][var_plscount] = v_edgeloc1
				w_decay_90[var_index][var_plscount] = v_edgelvl3
				w_decay_90_t[var_index][var_plscount] = v_edgeloc3
			endif
			edgestats/f=0.2/q/r=(var_peak_t,var_plscount*var_frq+var_end) $str_waves
			if (V_flag == 2) //failsafe
				variable/g var_fail = 1
			else //decay2080
				$("dc28_"+str_output)[var_index][var_plscount] = v_edgeloc3-v_edgeloc1
				w_decay_20[var_index][var_plscount] = v_edgelvl1
				w_decay_20_t[var_index][var_plscount] = v_edgeloc1
				w_decay_80[var_index][var_plscount] = v_edgelvl3
				w_decay_80_t[var_index][var_plscount] = v_edgeloc3
			endif
//AUC
			duplicate/o $str_waves, tempwave
			tempwave = tempwave-$("bs_"+str_output)[var_index][var_plscount]
			$("ar19_"+str_output)[var_index][var_plscount] = faverage(tempwave,var_r10,var_d90)
			killwaves tempwave
//width
			variable/g var_hm = 0.5*($("amp_"+str_output)[var_index][var_plscount])+$("bs_"+str_output)[var_index][var_plscount]
			findlevel/q/r=(var_peak_t,var_plscount*var_frq+var_start) $str_waves,var_hm
			variable/g var_hmrt = v_levelx
			w_fwhm1[var_index][var_plscount] = $str_waves(v_levelx)
			w_fwhm1_t[var_index][var_plscount] = v_levelx
			findlevel/q/r=(var_peak_t,var_plscount*var_frq+var_end) $str_waves,var_hm
			variable/g var_hmdt = v_levelx
			w_fwhm2[var_index][var_plscount] = $str_waves(v_levelx)
			w_fwhm2_t[var_index][var_plscount] = v_levelx
			$("hw_"+str_output)[var_index][var_plscount] = var_hmdt-var_hmrt
//multipeak
			if((numtype(str2num(str_multi)) == 0))
				wavestats/q/r=(var_plscount*var_frq+var_start-0.5*var_base_w,var_plscount*var_frq+var_start+0.5*var_base_w) $str_waves
				variable/g var_noise = v_sdev
				duplicate/o/r=(var_plscount*var_frq+var_start,var_plscount*var_frq+var_search) $str_waves w_work
				if (var_pol == 0)
					w_work = 0-w_work
				endif
				MPF2_AutoMPFit("temp_MP", "gauss", "peak_coefs_%d", "linear", "bsln", nameofwave(w_work), "", 5, noiseEst=var_noise, smFact=var_multi)
				killwaves/z W_AutoPeakInfo
				movewave root:temp_MP_0:W_AutoPeakInfo, root:
				if (dimsize(W_AutoPeakInfo,0)>0)
					make/o/n=(dimsize(W_AutoPeakInfo,0)) $("w_multi_pk_t_"+str_waves+"_"+num2str(var_plscount)) = nan
					$("w_multi_pk_t_"+str_waves+"_"+num2str(var_plscount)) = W_AutoPeakInfo[p][0]
					duplicate/o $("w_multi_pk_t_"+str_waves+"_"+num2str(var_plscount)) $("w_multi_pk_v_"+str_waves+"_"+num2str(var_plscount))
					$("w_multi_pk_v_"+str_waves+"_"+num2str(var_plscount)) = $str_waves($("w_multi_pk_t_"+str_waves+"_"+num2str(var_plscount)))
					appendtograph $("w_multi_pk_v_"+str_waves+"_"+num2str(var_plscount)) vs $("w_multi_pk_t_"+str_waves+"_"+num2str(var_plscount))
					modifygraph mode($("w_multi_pk_v_"+str_waves+"_"+num2str(var_plscount)))=8,marker($("w_multi_pk_v_"+str_waves+"_"+num2str(var_plscount)))=8,mrkThick($("w_multi_pk_v_"+str_waves+"_"+num2str(var_plscount)))=2,zcolor($("w_multi_pk_v_"+str_waves+"_"+num2str(var_plscount)))={w_colors[0,dimsize(W_AutoPeakInfo,0)],*,*,pastels,0}
					$("Npk_"+str_output)[var_index][var_plscount] = dimsize(W_AutoPeakInfo,0)
					$("multi_amp_"+str_output)[var_index][var_plscount][0,dimsize(W_AutoPeakInfo,0)-1] = $("w_multi_pk_v_"+str_waves+"_"+num2str(var_plscount))[r]-$("bs_"+str_output)[var_index][var_plscount]
					$("multi_pk_"+str_output)[var_index][var_plscount][0,dimsize(W_AutoPeakInfo,0)-1] = $("w_multi_pk_v_"+str_waves+"_"+num2str(var_plscount))[r]
					$("multi_pkT_"+str_output)[var_index][var_plscount][0,dimsize(W_AutoPeakInfo,0)-1] = $("w_multi_pk_t_"+str_waves+"_"+num2str(var_plscount))[r]-(var_plscount*var_frq+var_stim)
				else
					$("Npk_"+str_output)[var_index][var_plscount] = 0
				endif
				killwaves/z w_work
				killdatafolder/z $"temp_MP_0"
			endif
		else //spike
			$("APpb_"+str_output)[var_index][var_plscount] = 1
			$("APt_"+str_output)[var_index][var_plscount] = var_peak_t-(var_plscount*var_frq+var_stim)
			if (var_pol == 1) //upward
				w_peak_ap_t[var_index][var_plscount] = v_maxloc
				w_peak_ap[var_index][var_plscount] = v_max
			else //downward
				w_peak_ap_t[var_index][var_plscount] = v_minloc
				w_peak_ap[var_index][var_plscount] = v_min
			endif
		endif
//dV/dt
		if((numtype(str2num(str_dvdt)) == 0))
			if(var_plscount == 0)
				duplicate/o $str_waves $(str_waves+"_dvdt")
				differentiate $(str_waves+"_dvdt")
				appendtograph/r $(str_waves+"_dvdt")
			endif
			wavestats/q/r=(var_plscount*var_frq+var_start-0.5*var_base_w,var_plscount*var_frq+var_start+0.5*var_base_w) $(str_waves+"_dvdt")
			variable/g var_bs = v_avg
			$("dvdt_bs_"+str_output)[var_index][var_plscount] = v_avg
			$("dvdt_bsT_"+str_output)[var_index][var_plscount] = var_plscount*var_frq+var_start
			w_dvdt_base[var_index][var_plscount] = $(str_waves+"_dvdt")(var_plscount*var_frq+var_start)
			w_dvdt_base_t[var_index][var_plscount] = var_plscount*var_frq+var_start
			wavestats/q/r=(var_plscount*var_frq+var_start,var_plscount*var_frq+var_search) $(str_waves+"_dvdt")
			if (var_pol == 1) //upward
				variable/g var_peak_t = v_maxloc
				variable/g var_peak = v_max
			else //downward
				variable/g var_peak_t = v_minloc
				variable/g var_peak = v_min
			endif
			if (round(var_peak_t - (var_plscount*var_frq+var_start))<1) //failsafe
				var_peak_t = var_peak_t+1
				variable/g var_fail = 1
			endif
			if ((var_pol == 1 && var_peak-var_bs < var_dvdt) || (var_pol == 0 && var_peak-var_bs > var_dvdt)) //subthreshold
				wavestats/q/r=(var_peak_t-0.5*var_peak_w,var_peak_t+0.5*var_peak_w) $(str_waves+"_dvdt")
				$("dvdt_amp_"+str_output)[var_index][var_plscount] = v_avg-$("dvdt_bs_"+str_output)[var_index][var_plscount]
				$("dvdt_pk_"+str_output)[var_index][var_plscount] = v_avg
				$("dvdt_APpb_"+str_output)[var_index][var_plscount] = 0
				$("dvdt_pkT_"+str_output)[var_index][var_plscount] = var_peak_t-(var_plscount*var_frq+var_stim)
				if (var_pol == 1) //upward
					w_dvdt_peak_t[var_index][var_plscount] = v_maxloc
					w_dvdt_peak[var_index][var_plscount] = v_max
				else //downward
					w_dvdt_peak_t[var_index][var_plscount] = v_minloc
					w_dvdt_peak[var_index][var_plscount] = v_min
				endif
			else //spike
				$("dvdt_APpb_"+str_output)[var_index][var_plscount] = 1
				$("dvdt_APt_"+str_output)[var_index][var_plscount] = var_peak_t-(var_plscount*var_frq+var_stim)
				wavestats/q/r=(var_peak_t-0.5*var_peak_w,var_peak_t+0.5*var_peak_w) $(str_waves+"_dvdt")
				$("dvdt_APamp_"+str_output)[var_index][var_plscount] = v_avg-$("dvdt_bs_"+str_output)[var_index][var_plscount]
				if (var_pol == 1) //upward
					w_dvdt_peak_ap_t[var_index][var_plscount] = v_maxloc
					w_dvdt_peak_ap[var_index][var_plscount] = v_max
				else //downward
					w_dvdt_peak_ap_t[var_index][var_plscount] = v_minloc
					w_dvdt_peak_ap[var_index][var_plscount] = v_min
				endif
			endif
		endif
		variable/g var_scroll = 0
		do
			setdimlabel 0,var_index,$str_waves,$stringfromlist(var_scroll,str_output_list)
			var_scroll += 1
		while (var_scroll < itemsinlist(str_output_list))
		if(numtype(str2num(str_dvdt)) == 0)
			variable/g var_scroll = 0
			do
				setdimlabel 0,var_index,$str_waves,$stringfromlist(var_scroll,str_dvdt_list)
				var_scroll += 1
			while (var_scroll < itemsinlist(str_dvdt_list))
		endif
		if(numtype(str2num(str_multi)) == 0)
			variable/g var_scroll = 0
			do
				setdimlabel 0,var_index,$str_waves,$stringfromlist(var_scroll,str_multi_list)
				var_scroll += 1
			while (var_scroll < itemsinlist(str_multi_list))
		endif
		if (var_fail == 1) //failure
			variable/g var_scroll = 0
			do
				$stringfromlist(var_scroll,str_output_list)[var_index][var_plscount] = NaN
				var_scroll += 1
			while (var_scroll < itemsinlist(str_output_list))
			variable/g var_scroll = 0
			do
				$stringfromlist(var_scroll,str_display_list)[var_index][var_plscount] = NaN
				var_scroll += 1
			while (var_scroll < itemsinlist(str_display_list))
		endif
		variable/g var_scroll = 0
		do
			setdimlabel 1,var_plscount,$("pulse#"+num2str(var_plscount+1)),$stringfromlist(var_scroll,str_output_list)
			var_scroll += 1
		while (var_scroll < itemsinlist(str_output_list))
		if(numtype(str2num(str_dvdt)) == 0)
			variable/g var_scroll = 0
			do
				setdimlabel 1,var_plscount,$("pulse#"+num2str(var_plscount+1)),$stringfromlist(var_scroll,str_dvdt_list)
				var_scroll += 1
			while (var_scroll < itemsinlist(str_dvdt_list))
		endif
		if(numtype(str2num(str_multi)) == 0)
			variable/g var_scroll = 0
			do
				setdimlabel 1,var_plscount,$("pulse#"+num2str(var_plscount+1)),$stringfromlist(var_scroll,str_multi_list)
				var_scroll += 1
			while (var_scroll < itemsinlist(str_multi_list))
		endif
		var_plscount += 1
	while (var_plscount < var_pulse)
	var_index += 1
while (var_index < var_total)
//displaying landmarks
appendtograph w_base vs w_base_t
appendtograph w_peak vs w_peak_t
appendtograph w_onset vs w_onset_t
appendtograph w_rise_10 vs w_rise_10_t
appendtograph w_rise_90 vs w_rise_90_t
appendtograph w_rise_20 vs w_rise_20_t
appendtograph w_rise_80 vs w_rise_80_t
appendtograph w_decay_10 vs w_decay_10_t
appendtograph w_decay_90 vs w_decay_90_t
appendtograph w_decay_20 vs w_decay_20_t
appendtograph w_decay_80 vs w_decay_80_t
appendtograph w_fwhm1 vs w_fwhm1_t
appendtograph w_fwhm2 vs w_fwhm2_t
appendtograph w_peak_ap vs w_peak_ap_t
modifygraph mode(w_base)=3,mrkThick(w_base)=2,rgb(w_base)=(52428,52425,1)
modifygraph mode(w_peak)=3,mrkThick(w_peak)=2,rgb(w_peak)=(65535,0,0)
modifygraph mode(w_onset)=3,mrkThick(w_onset)=2,rgb(w_onset)=(65535,0,52428)
modifygraph mode(w_rise_10)=3,mrkThick(w_rise_10)=2,rgb(w_rise_10)=(2,39321,1)
modifygraph mode(w_rise_90)=3,mrkThick(w_rise_90)=2,rgb(w_rise_90)=(2,39321,1)
modifygraph mode(w_rise_20)=3,mrkThick(w_rise_20)=2,rgb(w_rise_20)=(16385,28398,65535)
modifygraph mode(w_rise_80)=3,mrkThick(w_rise_80)=2,rgb(w_rise_80)=(16385,28398,65535)
modifygraph mode(w_decay_10)=3,mrkThick(w_decay_10)=2,rgb(w_decay_10)=(52428,34958,1)
modifygraph mode(w_decay_90)=3,mrkThick(w_decay_90)=2,rgb(w_decay_90)=(52428,34958,1)
modifygraph mode(w_decay_20)=3,mrkThick(w_decay_20)=2,rgb(w_decay_20)=(36873,14755,58982)
modifygraph mode(w_decay_80)=3,mrkThick(w_decay_80)=2,rgb(w_decay_80)=(36873,14755,58982)
modifygraph mode(w_fwhm1)=3,mrkThick(w_fwhm1)=2,rgb(w_fwhm1)=(65535,65535,0)
modifygraph mode(w_fwhm2)=3,mrkThick(w_fwhm2)=2,rgb(w_fwhm2)=(65535,65535,0)
modifygraph mode(w_peak_ap)=3,mrkThick(w_peak_ap)=2,rgb(w_peak_ap)=(34952,34952,34952)
killwaves/z W_AutoPeakInfo
killwaves/z w_colors
//storing data
variable/g var_scroll = 0
do
	$("data_array_"+str_output)[][][var_scroll] = $stringfromlist(var_scroll,str_output_list)[p][q]
	copydimlabels/rows=0 $stringfromlist(var_scroll,str_output_list),$("data_array_"+str_output)
	copydimlabels/cols=1 $stringfromlist(var_scroll,str_output_list),$("data_array_"+str_output)
	setdimlabel 2,var_scroll,$stringfromlist(var_scroll,str_output_list),$("data_array_"+str_output)
	var_scroll += 1
while (var_scroll < itemsinlist(str_output_list))
if(numtype(str2num(str_multi)) == 0)
	variable/g var_scroll = 0
	do
		redimension/n=(dimsize($stringfromlist(var_scroll,str_multi_list),0),dimsize($stringfromlist(var_scroll,str_multi_list),1),wavemax($("Npk_"+str_output))) $stringfromlist(var_scroll,str_multi_list)
		copydimlabels/rows=0 $("data_array_"+str_output),$stringfromlist(var_scroll,str_multi_list)
		copydimlabels/cols=1 $("data_array_"+str_output),$stringfromlist(var_scroll,str_multi_list)
		variable/g var_label = 0
		do
			setdimlabel 2,var_label,$("peak#"+num2str(var_label+1)),$stringfromlist(var_scroll,str_multi_list)
			var_label+=1
		while(var_label<wavemax($("Npk_"+str_output)))
		var_scroll += 1
	while (var_scroll < itemsinlist(str_multi_list))
endif
if(numtype(str2num(str_dvdt)) == 0)
	appendtograph/r w_dvdt_base vs w_dvdt_base_t
	appendtograph/r w_dvdt_peak vs w_dvdt_peak_t
	appendtograph/r w_dvdt_peak_ap vs w_dvdt_peak_ap_t
	modifygraph mode(w_dvdt_base)=3,mrkThick(w_dvdt_base)=2,rgb(w_dvdt_base)=(34952,34952,34952)
	modifygraph mode(w_dvdt_peak)=3,mrkThick(w_dvdt_peak)=2,rgb(w_dvdt_peak)=(3,52428,1)
	modifygraph mode(w_dvdt_peak_ap)=3,mrkThick(w_dvdt_peak_ap)=2,rgb(w_dvdt_peak_ap)=(0,0,0)
	//variable/g var_r_min = min(wavemin(w_dvdt_base),wavemin(w_dvdt_peak),wavemin(w_dvdt_peak_ap))
	//variable/g var_r_max = max(wavemax(w_dvdt_base),wavemax(w_dvdt_peak),wavemax(w_dvdt_peak_ap))
	make/o/n=3 w_temp_min = {wavemin(w_dvdt_base),wavemin(w_dvdt_peak),wavemin(w_dvdt_peak_ap)}
	make/o/n=3 w_temp_max = {wavemax(w_dvdt_base),wavemax(w_dvdt_peak),wavemax(w_dvdt_peak_ap)}
	variable/g var_r_min = wavemin(w_temp_min)
	variable/g var_r_max = wavemax(w_temp_max)
	killwaves/z w_temp_min,w_temp_max
	if(var_r_min < 0)
		var_r_min = 1.1*var_r_min
	else
		var_r_min = 0.9*var_r_min
	endif
	if(var_r_max < 0)
		var_r_max = 0.9*var_r_max
	else
		var_r_max = 1.1*var_r_max
	endif
	var_r_max = max(abs(var_r_max),abs(var_r_min))
	var_r_min = min(0-abs(var_r_max),0-abs(var_r_min))
	setaxis right,var_r_min,var_r_max
	modifyGraph axisenab(right)={0.5,1}
	modifyGraph axisenab(left)={0,0.5}
	//variable/g var_binL = floor(min(wavemin($("dvdt_amp_"+str_output)),wavemin($("dvdt_APamp_"+str_output))))-1
	//variable/g var_binH = ceil(max(wavemax($("dvdt_amp_"+str_output)),wavemax($("dvdt_APamp_"+str_output))))+1
	make/o/n=2 w_temp_binL = {wavemin($("dvdt_amp_"+str_output)),wavemin($("dvdt_APamp_"+str_output))}
	make/o/n=2 w_temp_binH = {wavemax($("dvdt_amp_"+str_output)),wavemax($("dvdt_APamp_"+str_output))}
	variable/g var_binL = floor(wavemin(w_temp_binL))-1
	variable/g var_binH = ceil(wavemax(w_temp_binH))+1
	killwaves/z w_temp_binL,w_temp_binH
	//variable/g var_binN = 0.2*dimsize($("dvdt_amp_"+str_output),0)*dimsize($("dvdt_amp_"+str_output),1)
	variable/g var_binN = var_binH-var_binL
	make/o/n=(var_binN) $("dvdt_amp_"+str_output+"_histo")
	histogram/p/c/b={var_binL,1,var_binN} $("dvdt_amp_"+str_output), $("dvdt_amp_"+str_output+"_histo")
	make/o/n=(var_binN) $("dvdt_APamp_"+str_output+"_histo")
	histogram/p/c/b={var_binL,1,var_binN} $("dvdt_APamp_"+str_output), $("dvdt_APamp_"+str_output+"_histo")
	getwindow kwtopwin,wsize
	variable/g var_left_dvdt = v_right + 10
	variable/g var_top_dvdt = v_top
	variable/g var_right_dvdt = var_left_dvdt + 200
	variable/g var_bottom_dvdt = var_top_dvdt + 350
	display $("dvdt_amp_"+str_output+"_histo"),$("dvdt_APamp_"+str_output+"_histo")
	label bottom "dV/dt peak amplitude"
	label left "normalized counts"
	movewindow var_left_dvdt,var_top_dvdt,var_right_dvdt,var_bottom_dvdt
	modifygraph mode($("dvdt_amp_"+str_output+"_histo"))=7,hbFill($("dvdt_amp_"+str_output+"_histo"))=2,rgb($("dvdt_amp_"+str_output+"_histo"))=(0,0,0,32768)
	modifygraph mode($("dvdt_APamp_"+str_output+"_histo"))=7,hbFill($("dvdt_APamp_"+str_output+"_histo"))=2,rgb($("dvdt_APamp_"+str_output+"_histo"))=(65535,0,0,32768)
	legend/c/n=lgd/j/b=1/a=lt/f=0/x=0/y=0 "\\s(dvdt_amp_"+str_output+"_histo) subthreshold\r\\s(dvdt_APamp_"+str_output+"_histo) suprathreshold"
	setdrawenv xcoord=bottom,ycoord=left,linefgc=(34952,34952,34952),dash=3,linethick=2.00
	drawline var_dvdt,1,var_dvdt,0
endif
resumeupdate
end