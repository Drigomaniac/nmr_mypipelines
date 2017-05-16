%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Default Parameters; 
obj.Postproc.fx_cline.in.tmp = '';
obj.Postproc.trkland.fx.tmp.fx_roa = '';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%@ dwi_ADRC parameters (inits)
obj.Postproc.trkland.addpath = '<GET LOCATION OF TRK_LANDS>' ;
obj.Postproc.trkland.trks.fx.tmp.b0 = [filesep 'space' filesep 'public_html' ... 
    filesep 'rdp20' filesep 'fornix_ROA' filesep 'FX_1.8mm_orig' filesep ...
    '141219_8CS00178_b0.nii.gz' ] ;

obj.Postproc.trkland.trks.fx.tmp.fx_roa = [filesep 'space' filesep 'public_html' ... 
    filesep 'rdp20' filesep 'fornix_ROA' filesep 'FX_1.8mm_orig' filesep ...
    'ROA178noAC_141210_8CS00178_b0.nii.gz' ] ;

obj.Postproc.trkland.trks.fx.tmp.fx_tmp = [filesep 'space' filesep 'public_html' ... 
    filesep 'rdp20' filesep 'fornix_ROA' filesep 'FX_1.8mm_orig' filesep ...
    'TMP_178_wholeFornix.nii.gz' ] ;

function postproc_trkland (obj)
%Adding path if necessary:;
addpath(genpath(obj.Postproc.trkland.addpath));
trks_names = fields(obj.Postproc.trkland.trks);
for ii=1:numel(trks_names)
   %Init output variables:
   
end
end