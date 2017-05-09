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
        
        FS_location='/eris/bang/ADRC/FreeSurfer6.0/';
        
        
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
            
            
            %%%%%%%%%%%%
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
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For BET2:
            obj.Params.Bet2.in.movefiles = ['..' filesep '04_Bet'];
            obj.Params.Bet2.in.fracthrsh = 0.4;
            obj.Params.Bet2.in.fn = obj.Params.B0MoCo.out.fn;
            
            obj.proc_bet2();
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For EDDY:
            obj.Params.Eddy.in.movefiles = ['..' filesep '04_Eddy'];
            obj.Params.Eddy.in.fn=obj.Params.B0MoCo.out.fn;
            obj.Params.Eddy.in.bvals=obj.Params.B0MoCo.out.bvals';
            obj.Params.Eddy.in.bvecs=obj.Params.B0MoCo.out.bvecs;
            obj.Params.Eddy.in.mask = obj.Params.Bet2.out.mask;
            obj.Params.Eddy.in.index= ones(1,67); %for 35 volumes
            obj.Params.Eddy.in.acqp= [ 0 -1 0 0.08201 ]; %PE=A>>P (-1 at Y ) Echo spacing = 0.59 and EPI factor = 140 ==> 0.59(EPI)*0.001*(PE-1)139
            
            obj.proc_eddy();
            
            
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