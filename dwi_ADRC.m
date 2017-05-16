classdef dwi_ADRC < dwiMRI_Session
    %%  classdef dwi_ADRC < dwiMRI_Session
    %%  This class is a subclass of its parent class dwi_MRI_Session.m
    %%  (where it will inherent other methods).
    %%  Created by:
    %%              Aaron Schultz
    %%              Rodrigo D. Perea rperea@mgh.harvard.edu
    %%
    %%
    
    properties
        %root directoy where raw data lives:
        session_location='/eris/bang/ADRC/Sessions/';
        dcm_location = '/eris/bang/ADRC/DICOM_Archive/';
        gradfile='/autofs/space/kant_004/users/ConnectomeScanner/Scripts/adrc_diff_prep/bash/gradient_nonlin_unwarp/gradient_coil_files/coeff_AS302.grad';
        sh_gradfile=[ '/eris/bang/ADRC/Scripts/DEPENDENCIES/GradNonLin_Correc/run_mris_gradient_nonlin__unwarp_volume__batchmode_ADRC_v3.sh ' ...
            '/usr/pubsw/common/matlab/8.5 '];
      
        %sh dependencies:
        b0MoCo_rotate_bvecs_sh='/eris/bang/ADRC/Scripts/DEPENDENCIES/rotate_bvecs.sh'; %For rotating the bvecs after proc_b0MoCo
        init_rotate_bvecs_sh='/eris/bang/ADRC/Scripts/DEPENDENCIES/mod_fdt_rotate_bvecs.sh'; %For standarizing the bvecs after proc_dcm2nii
        col2rows_sh='/eris/bang/ADRC/Scripts/DEPENDENCIES/drigo_col2rows.sh';
        dependencies_dir='/eris/bang/ADRC/Scripts/DEPENDENCIES/';
        %FreeSurfer Dependencies
        FS_location='/eris/bang/ADRC/FreeSurfer6.0/';
        init_FS = '/usr/local/freesurfer/stable6';
        
    end
    methods
        function obj = dwi_ADRC(sessionname,opt)
             %For compiler code:
            if ~isdeployed()
                addpath(genpath('/autofs/space/kant_004/users/rdp20/scripts/matlab'));
            end
            
            %%%  If opt is passed, then the root Sessions folder will be
            %%%  replaced with this argument.
            if nargin>1
                obj.root = opt;
            end
            obj.sessionname = sessionname;
            obj.root = [obj.session_location sessionname '/DWIs/'];
            obj.dcm_location = [ obj.dcm_location sessionname filesep ];
            obj.session_location= [ obj.session_location sessionname filesep ] ;
           
            %If the folder /DWIs/ does not exist, then create it!
            if exist(obj.root,'dir')==0
                obj.make_root();
            end
            obj.objectHome = obj.root ; 
            %Check to see if a *.mat file exists.
            if exist([obj.objectHome filesep sessionname '.mat'],'file')>0
                load([obj.objectHome filesep sessionname '.mat']);
                oldroot = obj.root;
                obj.wasLoaded = true;
            else
                obj.setMyParams;
            end
            
            %Check if *.nii.gz files exist, if not get them from DCM2nii:
            obj.rawfiles = dir_wfp([obj.root 'Orig/*.nii.gz' ] );
            if isempty(obj.rawfiles) || numel(obj.rawfiles) ~= 4 % 4 DWIs sequence acquired here
                RunFlag = true ;
                obj.getDCM2nii(RunFlag);
            end
            
            if nargin>1
                if ~strcmpi(oldroot,newroot)
                    obj = replaceObjText(obj,{oldroot},{newroot});
                    obj.resave;
                end
            end
            %Start the CommonProc:;;
            obj.CommonProc();
        end
        
        function obj=setMyParams(obj)
            %%%%%%%%%%%%
            %Global parameters:
            obj.vox = [1.8 1.8 1.8];
            obj.setDefaultParams; %this will call the method in the superclass dwiMRI_Session.m 
            obj.rawfiles = dir_wfp([obj.root 'Orig' filesep '*.nii.gz' ] );
            
        end
        
        function obj = CommonProc(obj)
            obj.dosave = true ; %To record process in MAT file
            
             if isempty(obj.rawfiles)
                obj.rawfiles = dir_wfp([obj.root 'Orig_' filesep '*.nii' ] );
            end
            
            %%%%%%%%%%%%
            %For proc_dropvols
            obj.Params.DropVols.in.tmin='1';
            obj.Params.DropVols.in.tsize='67';
            obj.Params.DropVols.in.prefix='dv_';
            obj.Params.DropVols.in.movefiles=['..' filesep '01_DropVols' filesep ];
            obj.Params.DropVols.in.fn=obj.rawfiles;
            %Bvecs and bvals will be created from XX.in.fn and XX.out.fn
            obj.Params.DropVols.out.fn=dir_wfp([obj.root, '01_DropVols', filesep, '*.nii.gz']);
            
            obj.proc_drop_vols();
            
            %%%%%%%%%%%%
            %For gradient non-linearity correction
            obj.Params.GradNonlinCorrect.in.movefiles = '../02_GradCorrect/';
            obj.Params.GradNonlinCorrect.in.prefix = 'gnc_';
            obj.Params.GradNonlinCorrect.in.gradfile = obj.gradfile;
            
            obj.Params.GradNonlinCorrect.in.fslroi = [ 0 1 ]; %To extraact the 1st b0
            obj.Params.GradNonlinCorrect.in.fn = obj.Params.DropVols.out.fn;

            obj.proc_gradient_nonlin_correct();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For FreeSurfer Segmentation (needed for b0 bbreg correction)
             obj.Params.FreeSurfer.dir = [ filesep 'eris' filesep 'bang' filesep ...
            'ADRC' filesep 'FreeSurfer6.0' ] ;
            obj.Params.FreeSurfer.init_location = obj.init_FS;
            
        
            %Retrieving a T1 scan:
            [sys_error, obj.Params.FreeSurfer.in.T1 ] = system(['ls ' obj.session_location 'T1' filesep '*.nii | head -1' ]);
            if sys_error ~= 0 %No problem, we get the T1 the continue...
                fprintf(['\nError when finding the T1:'  obj.Params.FreeSurfer.in.T1  '\n'])
            end
            
            %Retrieving a T2 scan:
            [sys_error, obj.Params.FreeSurfer.in.T2 ] = system(['ls ' obj.session_location 'other' filesep '*T2SPACE* | head -1' ]);
            if sys_error ~= 0 %No problem, we get the T1 the continue...
                fprintf(['\nNo T2 found:'  obj.Params.FreeSurfer.in.T2  '\n'])
                obj.Params.FreeSurfer.in.T2exist=false;
            else
                obj.Params.FreeSurfer.in.T2exist=true;
            end
            
            obj.Params.FreeSurfer.out.aparcaseg = [ obj.Params.FreeSurfer.dir ...
                filesep obj.sessionname filesep 'mri' filesep 'aparc+aseg.mgz' ] ;
            
            obj.proc_getFreeSurfer();
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For b0 motion correciton (based on interspersed b0s)
            obj.fsdir=[obj.FS_location obj.sessionname ] ; 
            obj.Params.B0MoCo.FS = obj.fsdir;
            obj.Params.B0MoCo.in.movefiles = '../03_B0s_MoCo/';
            obj.Params.B0MoCo.in.prefix = 'moco_';
            obj.Params.B0MoCo.in.nDoF = '12' ; 
            obj.Params.B0MoCo.in.grad_rel = obj.Params.GradNonlinCorrect.out.warpfile;
            
            obj.Params.B0MoCo.in.fn = obj.Params.GradNonlinCorrect.out.fn;
            obj.Params.B0MoCo.in.bvals = obj.Params.DropVols.out.bvals;
            obj.Params.B0MoCo.in.bvecs = obj.Params.DropVols.out.bvecs;
            obj.Params.B0MoCo.in.sh_rotate_bvecs = obj.b0MoCo_rotate_bvecs_sh; 
            obj.proc_b0s_MoCo();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For BET2:
            obj.Params.Bet2.in.movefiles = ['..' filesep '04_Bet'];
            obj.Params.Bet2.in.fracthrsh = 0.4;
            obj.Params.Bet2.in.fn = obj.Params.B0MoCo.out.fn;
            
            obj.proc_bet2();
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For EDDY:
            obj.Params.Eddy.in.movefiles = ['..' filesep '04_Eddy'];
            obj.Params.Eddy.in.fn=obj.Params.B0MoCo.out.fn;
            obj.Params.Eddy.in.bvals=obj.Params.B0MoCo.out.bvals';
            obj.Params.Eddy.in.bvecs=obj.Params.B0MoCo.out.bvecs;
            obj.Params.Eddy.in.mask = obj.Params.Bet2.out.mask;
            obj.Params.Eddy.in.index= ones(1,67); %for 35 volumes
            obj.Params.Eddy.in.acqp= [ 0 -1 0 0.08201 ]; %PE=A>>P (-1 at Y ) Echo spacing = 0.59 and EPI factor = 140 ==> 0.59(EPI)*0.001*(PE-1)139
            
            obj.proc_eddy();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %To generate a mask after Eddy:
            %This step will 1) define a better mask if eddy affected the
            %movement of the head and 2) remove issues known to happen at
            %the edges of the brain when using the --wls option in dtifit!\
            %Reference: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;6eb4d787.1610
            obj.Params.MaskAfterEddy.in.movefiles = ['..' filesep '05_MaskAfterEddy'];
            obj.Params.MaskAfterEddy.in.fn = obj.Params.Eddy.in.fn; %Since we don't have a b0, we pass the full dwi and the method will take care of
            obj.Params.MaskAfterEddy.in.prefix = 'after_eddy';
            obj.Params.MaskAfterEddy.in.fracthrsh = '0.4';
            
            obj.proc_mask_after_eddy();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Until now, we have treated all images separately, so now we
            %will rigid-body coregister all *using flirst and spline
            %interpolation (proved to be more accurate):
            obj.Params.CoRegMultiple.in.fn = obj.Params.Eddy.in.fn;
            obj.Params.CoRegMultiple.in.b0 = obj.Params.MaskAfterEddy.in.b0 ;
            obj.Params.CoRegMultiple.in.bvals = obj.Params.Eddy.in.bvals;
            obj.Params.CoRegMultiple.in.bvecs = obj.Params.Eddy.out.bvecs;
            
            obj.Params.CoRegMultiple.in.movefiles = ['..' filesep '06_CoRegDWIs'];
            obj.Params.CoRegMultiple.in.ref_iteration = 2; % All images will be registered to this iteration (in ADRC, 7p5_set1, index 1 is for 2p7_set4!)
            
            
            obj.proc_coreg_multiple();
            
            
            
            
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             %For DTIFIT:
%             %We will use the --wls option as it seems to improve the fit of
%             %the diffusion tensor model and negativity values due to noise
%             %REF: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;735f6320.1409
%             obj.Params.Dtifit.in.movefiles = [ '..' filesep '05_Dtifit' ];
%             obj.Params.Dtifit.in.fn = obj.Params. Eddy.out.fn;
%             obj.Params.Dtifit.in.prefix = 'DTIFIT_FSLv509' ; %Double check this so you prefix the version of FSL!
%             obj.Params.Dtifit.in.bvecs = obj.Params.Eddy.out.bvecs;
%             obj.Params.Dtifit.in.bvals = obj.Params.Eddy.in.bvals;
%             obj.Params.Dtifit.in.mask = obj.Params.MaskAfterEddy.out.finalmask;
%             
%             obj.proc_dtifit();
%             
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             %For GQI:
%             obj.Params.GQI.in.movefiles = [ '..' filesep '05_Recon_gqi' ];
%             obj.Params.GQI.in.fn = obj.Params. Eddy.out.fn;
%             obj.Params.GQI.in.bvecs = obj.Params.Eddy.out.bvecs;
%             obj.Params.GQI.in.bvals = obj.Params.Eddy.in.bvals;
%             obj.Params.GQI.in.mask = obj.Params.MaskAfterEddy.out.finalmask;
%             obj.Params.GQI.in.prefix = 'GQI_DSISv041917' ; %Double check this so you prefix the version of DSISTUDIO!
%             obj.Params.GQI.out.export = 'gfa,nqa0,nqa1';
%             
%             obj.proc_gqi();
%             
%             
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             %FS2dwi:
%             obj.Params.FS2dwi.in.movefiles = ['..' filesep '05_FS2dwi' ];
%             obj.Params.FS2dwi.in.b0 = obj.Params.B0mean.out.fn ; 
%             obj.Params.FS2dwi.in.aparcaseg = obj.Params.FreeSurfer.out.aparcaseg ; 
%             
%             obj.Params.FS2dwi.in.tmpfile_aparcaseg = [ obj.dependencies_dir 'FS_aparc.txt' ] ; 
%             obj.Params.FS2dwi.in.tmpfile_aparcaseg2009 = [ obj.dependencies_dir 'FS_aparc2009.txt' ] ; 
%             obj.Params.FS2dwi.in.tmpfile_hippo_bil = [ obj.dependencies_dir 'FS_hippolabels_bil.txt' ] ;
%             
%             
%             
%             obj.Params.FS2dwi.in.aparcaseg2009 = ...
%                 strtrim(strrep(obj.Params.FreeSurfer.out.aparcaseg,'aparc+aseg','aparc.a2009s+aseg')); 
%             
%             %A possible error is the naming convention when only a T1 was
%             %used!!
%             obj.Params.FS2dwi.in.hippofield_left = ...
%                 strtrim(strrep(obj.Params.FreeSurfer.out.aparcaseg,'aparc+aseg','lh.hippoSfLabels-T1-T2.v10.FSvoxelSpace')); 
%             obj.Params.FS2dwi.in.hippofield_right = ...
%                 strtrim(strrep(obj.Params.FreeSurfer.out.aparcaseg,'aparc+aseg','rh.hippoSfLabels-T1-T2.v10.FSvoxelSpace')); 
%           
%             obj.proc_FS2dwi();
            
            
        end
        
        function resave(obj)
            save([obj.objectHome filesep obj.sessionname '.mat'],'obj');
        end
    end
    
    methods ( Access = protected )
        function obj = getDCM2nii(obj,torun)
            %For proc_DCM2NII:
            obj.Params.DCM2NII.specific_vols=68;
            obj.Params.DCM2NII.scanlog = [ obj.session_location filesep 'LogFiles' ...
                filesep 'scan.log' ] ;
            if ~exist(obj.Params.DCM2NII.scanlog,'file')
                error(['No scanlog file found in: ' obj.Params.DCM2NII.scanlog ' . Exiting...']);
            end
            obj.Params.DCM2NII.seq_names={ 'ep2d_diff_7p5k_set1E60' 'ep2d_diff_7p5k_set2E60' ...
                'ep2d_diff_7p5k_set3E60' 'ep2d_diff_2p5k_set4E60' };
            
           
            for ii=1:4 % 4 sets of DWIs in this project!
                  obj.Params.DCM2NII.in(ii).fsl2std_param = '-1 0 0 250.199 \n0 1 0 250.199 \n0 0 -1 0 \n0 0 0 1';
           
                obj.Params.DCM2NII.in(ii).prefix = obj.Params.DCM2NII.seq_names(ii);
                [ ~ , obj.Params.DCM2NII.in(ii).nvols ] = system([ 'cat ' ...
                    obj.Params.DCM2NII.scanlog ' | grep ' obj.Params.DCM2NII.seq_names{ii} ...
                    ' | tail -1 | awk ''{ print $7 }'' ' ]);
                obj.Params.DCM2NII.in(ii).nvols=str2num(obj.Params.DCM2NII.in(ii).nvols);
                [ ~ , obj.Params.DCM2NII.in(ii).first_dcmfiles ]= system([ 'cat ' ...
                    obj.Params.DCM2NII.scanlog ' | grep ' obj.Params.DCM2NII.seq_names{ii} ...
                    ' | tail -1 | awk ''{ print $8 }'' ' ]);
                
                obj.Params.DCM2NII.out(ii).location = [ obj.root 'Orig' filesep ];
                obj.Params.DCM2NII.out(ii).fn = [  obj.Params.DCM2NII.out(ii).location  cell2char(obj.Params.DCM2NII.seq_names(ii)) '.nii.gz' ];
            end
            if (torun) ; obj.proc_dcm2nii ; end 
        end
    end
end