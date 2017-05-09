classdef dwi_HAB < dwiMRI_Session
    %%  classdef dwi_ADRC < dwiMRI_Session
    %%  This class is a subclass of its parent class dwi_MRI_Session.m
    %%  (where it will inherent other methods).
    %%  Created by:
    %%              Aaron Schultz aschultz@martinos.org
    %%              Rodrigo Perea rpereacamargo@mgh.harvard.edu
    %%
    %%
    %%      Dependencies:
    %%          -FreeSurfer v6.0
    %%          -SPM8
    %%          -Ants tools
    %%          -DSI_studio
    %%  *Only filesep wit '/' are used in the properties class declaration.
    %%   Besides these, all should be any operating system compatible (tested in CentOS Linux)
    properties
        
        %root directoy where raw data lives:
        root_location='/cluster/sperling/HAB/Project1/DWIs_30b700/Sessions/';
        dcm_location='/cluster/sperling/HAB/Project1/DICOM_ARCHIVE/All_DICOMS/';
        session_location='/cluster/sperling/HAB/Project1/Sessions/';
        
        %template dependencies:
        HABn272_meanFA='/cluster/hab/HAB/Project1/DWIs_30b700/DEPENDENCIES/HABn272_MNI_Target/HABn272_meanFA.nii.gz';
        HABn272_meanFA_skel_dst='/cluster/hab/HAB/Project1/DWIs_30b700/DEPENDENCIES/HABn272_MNI_Target/HABn272_meanFA_skeleton_mask_dst.nii.gz';
        ref_region='/usr/pubsw/packages/fsl/5.0.9/data/standard/LowerCingulum_1mm.nii.gz';
        
        %sh dependencies:
        init_rotate_bvecs_sh='/cluster/sperling/HAB/Project1/Scripts/DWIs/mod_fdt_rotate_bvecs.sh'; %For standarizing the bvecs after proc_dcm2nii
        
        %skel_TOI dependencies
        skeltoi_location='/cluster/hab/HAB/Project1/DWIs_30b700/DEPENDENCIES/edJHU_MASK_ROIs/';
        skeltoi_tois={ 'not_ATR_edJHU' 'not_CCG_edJHU' 'not_CHIPP_edJHU' 'not_CST_edJHU'  ...
            'not_Fma_edJHU' 'not_Fmi_edJHU' 'not_IFOF_edJHU' 'not_ILF_edJHU' ...
            'not_SLF_no_temp_edJHU'  'ATR_L_edJHU' 'ATR_R_edJHU' ...
            'CCG_L_edJHU' 'CCG_R_edJHU' 'CHIPP_L_edJHU' 'CHIPP_R_edJHU' 'CST_L_edJHU' ...
            'CST_R_edJHU' 'Fma_edJHU' 'Fmi_edJHU' 'IFOF_L_edJHU' 'IFOF_R_edJHU' ...
            'ILF_L_edJHU' 'ILF_R_edJHU' 'SLF_L_edJHU' 'SLF_L_edJHU_no_temporal' ...
            'SLF_R_edJHU' 'SLF_R_edJHU_no_temporal' 'allTracts_edJHU' 'Global_skel'  };
    end
    
    methods
        function obj = dwi_HAB(sessionname,opt)
            
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
            obj.root = [obj.root_location sessionname '/DWIs/'];
            obj.dcm_location= [ obj.dcm_location sessionname filesep ];
            obj.session_location= [ obj.session_location sessionname filesep ];
            
            %If the folder <XX>/DWIs/ does not exist, then create it!
            if exist(obj.root,'dir')==0
                obj.make_root();
            end
            obj.objectHome = obj.root ;
            if exist([obj.objectHome filesep sessionname '.mat'],'file') > 0
                load([obj.objectHome filesep sessionname '.mat']);
                oldroot = obj.root;
                obj.wasLoaded = true;
            else
                obj.setMyParams;
            end
            
            %Check if *.nii.gz files exist, if not get them from DCM2nii:
            obj.rawfiles = dir_wfp([obj.root 'Orig/*.nii.gz' ] );
            if isempty(obj.rawfiles)
                RunFlag=true;
                obj.getDCM2nii(RunFlag);
            end
            
            if nargin>1
                if ~strcmpi(oldroot,newroot)
                    obj = replaceObjText(obj,{oldroot},{newroot});
                    obj.resave;
                end
            end
            %Continue with CommonProc
            obj.CommonProc();
        end
        
        function obj=setMyParams(obj)
            %Global parameters:
            obj.vox= [2 2 2 ];
            obj.setDefaultParams; %from dwiMRI_Session class
        end
        
        function obj = CommonProc(obj)
            obj.dosave = true ; %To record process in MAT file
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Get DCM2NII File location:
            %
            RunFlag=false; %here, we only allocate variable, not run proc_dcm2nii
            obj.getDCM2nii(RunFlag);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For BET2:
            obj.Params.Bet2.in.movefiles = ['..' filesep '01_Bet'];
            obj.Params.Bet2.in.fracthrsh = 0.4;
            obj.Params.Bet2.in.fn = obj.rawfiles;
            
            obj.proc_bet2();
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For EDDY:
            obj.Params.Eddy.in.movefiles = ['..' filesep '02_Eddy'];
            obj.Params.Eddy.in.fn=obj.rawfiles;
            obj.Params.Eddy.in.bvals=strrep(obj.rawfiles,'.nii.gz','.bvals');
            obj.Params.Eddy.in.bvecs=strrep(obj.rawfiles,'.nii.gz','.voxel_space.bvecs');
            obj.Params.Eddy.in.mask = obj.Params.Bet2.out.mask;
            obj.Params.Eddy.in.index= ones(1,35) ; %for 35 volumes
            obj.Params.Eddy.in.acqp= [ 0 -1 0 0.102]; %based on HAB diff sequence
            
            obj.proc_eddy();
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %To calculate mean motion based on edddy output:
            obj.Params.EddyMotion.in.movefiles = ['..' filesep '03_EddyMotion']; 
            obj.Params.EddyMotion.in.fn_eddy = obj.Params.Eddy.out.fn ;
            
            obj.proc_get_eddymotion();
                        
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For B0mean:
            obj.Params.B0mean.in.movefiles = ['..' filesep '03_B0mean'];
            obj.Params.B0mean.in.fn=obj.Params.Eddy.out.fn;
            obj.Params.B0mean.in.b0_nvols=5;
            
            obj.proc_meanb0();
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %To generate a mask after Eddy:
            %This step will 1) define a better mask if eddy affected the
            %movement of the head and 2) remove issues known to happen at
            %the edges of the brain when using the --wls option in dtifit!\
            %Reference: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;6eb4d787.1610
            obj.Params.MaskAfterEddy.in.movefiles = ['..' filesep '04_MaskAfterEddy'];
            obj.Params.MaskAfterEddy.in.fn=obj.Params.B0mean.out.fn;
            obj.Params.MaskAfterEddy.in.prefix = 'after_eddy';
            obj.Params.MaskAfterEddy.in.fracthrsh = '0.4';
            
            obj.proc_mask_after_eddy();
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For DTIFIT:
            %We will use the --wls option as it seems to improve the fit of
            %the diffusion tensor model and negativity values due to noise
            %REF: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;735f6320.1409
            obj.Params.Dtifit.in.movefiles = [ '..' filesep '05_Recon_dtifit' ];
            obj.Params.Dtifit.in.fn = obj.Params. Eddy.out.fn;
            obj.Params.Dtifit.in.prefix = 'DTIFIT_FSLv509' ; %Double check this so you prefix the version of FSL!
            obj.Params.Dtifit.in.bvecs = obj.Params.Eddy.out.bvecs;
            obj.Params.Dtifit.in.bvals = obj.Params.Eddy.in.bvals;
            obj.Params.Dtifit.in.mask = obj.Params.MaskAfterEddy.out.finalmask;
            
            obj.proc_dtifit();
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For GQI:
            obj.Params.GQI.in.movefiles = [ '..' filesep '05_Recon_gqi' ];
            obj.Params.GQI.in.fn = obj.Params. Eddy.out.fn;
            obj.Params.GQI.in.bvecs = obj.Params.Eddy.out.bvecs;
            obj.Params.GQI.in.bvals = obj.Params.Eddy.in.bvals;
            obj.Params.GQI.in.mask = obj.Params.MaskAfterEddy.out.finalmask;
            obj.Params.GQI.in.prefix = 'GQI_DSISv041917' ; %Double check this so you prefix the version of DSISTUDIO!
            obj.Params.GQI.out.export = 'gfa,nqa0,nqa1';
            
            obj.proc_gqi();
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For AntsReg:
            obj.Params.AntsReg.in.movefiles = ['..' filesep '06_Ants_CoReg' ];
            obj.Params.AntsReg.in.fn = obj.Params.Dtifit.out.FA ;
            obj.Params.AntsReg.in.ref = obj.HABn272_meanFA;
            obj.Params.AntsReg.in.threads = '4' ;
            obj.Params.AntsReg.in.prefix = 'Antsv201_2_HABn272_' ;
            
            obj.proc_antsreg();
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For Skeletonize:
            obj.Params.Skeletonize.in.movefiles = ['..' filesep '07_Skeletonize' ];
            obj.Params.Skeletonize.in.fn = obj.Params.AntsReg.out.fn ;
            obj.Params.Skeletonize.in.meanFA = obj.HABn272_meanFA;
            obj.Params.Skeletonize.in.skel_dst = obj.HABn272_meanFA_skel_dst;
            obj.Params.Skeletonize.in.thr = '0.3';
            obj.Params.Skeletonize.in.ref_region = obj.ref_region;
            obj.Params.Skeletonize.in.prefix = [ 'FSLv509_skelHABn272' ] ; %check fsl versio you are calling!
            
            obj.proc_skeletonize();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For Skel_TOI:
            obj.Params.Skel_TOI.in.location = obj.skeltoi_location;
            obj.Params.Skel_TOI.in.masks = obj.skeltoi_tois ;
            obj.Params.Skel_TOI.in.ext = '.nii.gz' ;
            obj.Params.Skel_TOI.in.suffix = '_n272TMP';
           
            obj.proc_getskeltois();
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~isdeployed() %Not compiled code, so run this!
                %Uploading Skel Data into DataCentral:
                UploadData_DWI(obj)
                
                if obj.dosave
                    save([obj.objectHome filesep obj.sessionname '.mat'],'obj');
                end
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For FreeSurfer 
            obj.Params.FreeSurfer.dir = [ filesep 'eris' filesep 'bang' filesep ...
            'HAB_Project1' filesep 'FreeSurferv6.0' ] ;
            %Retrieving a T1 scan:
            [sys_error, obj.Params.FreeSurfer.in.T1 ] = system(['ls ' obj.session_location 'MPRAGE' filesep '*.mgz | head -1' ]);
            if sys_error ~= 0 %No problem, we get the T1 the continue...
                fprintf(['\nError when finding the T1:'  obj.Params.FreeSurfer.in.T1  '\n'])
            end
            
            %Retrieving a T2 scan:
            [sys_error, obj.Params.FreeSurfer.in.T2 ] = system(['ls ' obj.session_location 'other' filesep '*T2* | head -1' ]);
            if sys_error ~= 0 %No problem, we get the T1 the continue...
                fprintf(['\nNo T2 found:'  obj.Params.FreeSurfer.in.T2  '\n'])
                obj.Params.FreeSurfer.in.T2exist=false;
            else
                obj.Params.FreeSurfer.in.T2exist=true;
            end
            
            obj.Params.FreeSurfer.out.aparcaseg = [ obj.Params.FreeSurfer.dir ...
                filesep obj.sessionname filesep 'mri' filesep 'aparc+aseg.mgz' ] 
            
            obj.proc_getFreeSurfer();
        end
        
        function resave(obj)
            save([obj.objectHome filesep obj.sessionname '.mat'],'obj');
        end
    end
    
    methods ( Access = protected )
        function obj = getDCM2nii(obj,torun)
            %For proc_DCM2NII:
            obj.Params.DCM2NII.specific_vols=35;
            obj.Params.DCM2NII.scanlog = [ obj.session_location  filesep 'LogFiles' ...
                filesep 'scan.log' ] ;
            
            if ~exist(obj.Params.DCM2NII.scanlog,'file')
                error(['No scanlog found in:' obj.Params.DCM2NII.scanlog '. Exiting...']);
            end
            
            obj.Params.DCM2NII.seq_names='DIFFUSION_HighRes_30';
            try
                %we modified grep err instead of tail -1 due to error in some subjects (e.g. 120419_4TT01420)
                exec_cmd=[ 'cat ' obj.Params.DCM2NII.scanlog ' | grep ' obj.Params.DCM2NII.seq_names ' | grep " 35 " | grep err | awk ''{ print $7 }'' ' ];
                [ ~ , obj.Params.DCM2NII.in.nvols ] = system(exec_cmd);
            catch
                errormsg=['DCM2NII: No 35 vols. when reading scanlog located in: ' ...
                    obj.Params.DCM2NII.scanlog '\n' ];
                obj.UpdateErrors(errormsg);
            end
            obj.Params.DCM2NII.in.fsl2std_param = '-1 0 0 254 \n0 1 0 254 \n0 0 -1 0 \n0 0 0 1';
            for ii=1:1 %For object compatiblity with ADRC that contains 3 DWIs sequences
                obj.Params.DCM2NII.in(ii).nvols=str2double(obj.Params.DCM2NII.in.nvols);
                %we modified grep err instead of tail -1 due to error in some subjects (e.g. 120419_4TT01420)
                exec_cmd=[ 'cat ' obj.Params.DCM2NII.scanlog ' | grep ' obj.Params.DCM2NII.seq_names ' | grep " 35 " | grep err | awk ''{ print $8 }'' ' ];
                [ ~ , obj.Params.DCM2NII.in(ii).first_dcmfiles ] = system(exec_cmd);
                
                obj.Params.DCM2NII.out(ii).location = [ obj.root 'Orig' filesep ];
                obj.Params.DCM2NII.out(ii).fn = [ obj.Params.DCM2NII.out(ii).location obj.Params.DCM2NII.seq_names '.nii.gz' ];
                
            end
            if (torun) ; obj.proc_dcm2nii ; end
        end
    end
end