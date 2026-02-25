#pragma rtGlobals=1

Macro Get_Waves_In_Folders_KO()

pGWIFko()

proc pGWIFko(str_input)
string str_input
prompt str_input "input waves ?"
variable/g var_Fcount = countobjects("",4)
string/g str_Flist
if(itemsinlist(str_Flist)!=0)
	killstrings str_Flist
endif
string/g str_Flist
variable/g var_Fscroll = 0
silent 1
pauseupdate
do
	//str_Flist = str_Flist+";"+getindexedobjname("", 4, var_Fscroll)
	str_Flist = addlistitem(getindexedobjname("", 4, var_Fscroll),str_Flist,";",var_Fscroll)
	var_Fscroll+=1
while(var_Fscroll<var_Fcount)
variable/g var_Fscroll = 0
do
	string/g str_Fbrowse = stringfromlist(var_Fscroll,str_Flist)
	//if(strlen(str_Fbrowse)==0)
	//	break
	//	var_Fscroll +=1
	//endif
	setdatafolder(str_Fbrowse)
	string/g str_wavelist = wavelist(str_input,";","")
	variable/g var_Wscroll = 0
	do
		if(strlen(str_wavelist)==0)
			break
		endif
		//string str_wave = getindexedobjname(str_Fbrowse,1,var_Wscroll)
		string/g str_wave = stringfromlist(var_Wscroll,str_wavelist)
		duplicate/o $str_wave, $(nameofwave($str_wave)+"_d")
		movewave $(nameofwave($str_wave)+"_d"), root:$(str_wave+getdatafolder(0)[1,strlen(getdatafolder(0))-2])
		var_Wscroll+=1
	while(var_Wscroll<itemsinlist(str_wavelist))
	setdatafolder root:
	var_Fscroll+=1
while(var_Fscroll<var_Fcount)
end