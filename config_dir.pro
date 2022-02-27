; 
; Initialize useful directories.
; 
; James Ruppert
; 10/26/18
; 
pro config_dir, dirs=dirs

scdir='/scratch/06040/tg853394/'
wkdir='/work2/06040/tg853394/stampede2/'
home='/home1/06040/tg853394/'
figdir='/home1/06040/tg853394/idl/figures/'

ctfil='/home1/06040/tg853394/idl/code/git/idlcode/misc/ctbl_jhr_17Apr2019.tbl'

dirs={scdir:scdir,wkdir:wkdir,home:home,figdir:figdir,ctfil:ctfil}

end
