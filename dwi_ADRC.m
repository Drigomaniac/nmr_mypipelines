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
       root = '/eris/bang/ADRC/Sessions/';
       dcm_location = '/eris/bang/ADRC/DICOM_Archive/';
       session_location='/eris/bang/ADRC/Sessions/';
       gradfile='/autofs/space/kant_004/users/ConnectomeScanner/Scripts/adrc_diff_prep/bash/gradient_nonlin_unwarp/gradient_coil_files/coeff_AS302.grad';
       sh_gradfile=[ '/autofs/space/kant_004/users/ConnectomeScanner/Scripts/adrc_diff_prep/run_mris_gradient_nonlin__unwarp_volume__batchmode_ADRC_v3.sh ' ...
                    '/usr/pubsw/common/matlab/8.5 '];
    end
    methods
        function obj = dwi_ADRC(sessionname,opt)
            %%%  If opt is passed, then the root Sessions folder will be 
            %%%  replaced with this argument.
            if nargin>1
                obj.root = opt;
            end
            
            obj.sessionname = sessionname;
            obj.root = [obj.root sessionname '/DWIs/'];
            obj.dcm_location = [ obj.dcm_location sessionname filesep ];
            obj.session_location= [ obj.session_location sessionname filesep  ];
            %If the folder /DWIs/ does not exist, then create it!             
            if exist(obj.root,'dir')==0
                try
                    system(['mkdir ' obj.root ]); 
                catch
                    disp([ 'Trying to create /DWIs/ in:' obj.root ... 
                        ' .Does the parent folder exist?'])
                end
            end
            obj.objectHome = obj.root;
            
            %newroot = obj.root;
            %oldroot = obj.root;
            %obj.setSPM12;  %No needed yet
            %obj.dosave = true;

            if exist([obj.objectHome filesep sessionname '.mat'],'file')>0 
                load([obj.objectHome filesep sessionname '.mat']);
                oldroot = obj.root;
                obj.wasLoaded = true;
            else
                obj.setMyParams; 
            end
            
            if nargin>1
                if ~strcmpi(oldroot,newroot)
                    obj = replaceObjText(obj,{oldroot},{newroot});
                    obj.resave;
                end
            end
        end
        
        function obj=setMyParams(obj)
            %%%%%%%%%%%%
            %Global parameters:
            obj.vox = [1.8 1.8 1.8];
            obj.rawfiles = dir_wfp([obj.root 'Orig' filesep '*.nii.gz' ] );
            if isempty(obj.rawfiles)
                obj.rawfiles = dir_wfp([obj.root 'Orig' filesep '*.nii' ] );
            end
            
            %%%%%%%%%%%%
            %For proc_DCM2NII:
            obj.Params.DCM2NII.specific_vols=68;
            obj.Params.DCM2NII.scanlog = [ obj.session_location filesep 'LogFiles' ...
                filesep 'scan.log' ] ;
            if ~exist(obj.Params.DCM2NII.scanlog,'file')
                error(['No scanlog file found in: ' obj.Params.DCM2NII.scanlog ' . Exiting...']);
            end
            dwi_setnames={ 'ep2d_diff_7p5k_set1E60' 'ep2d_diff_7p5k_set2E60' ...
                'ep2d_diff_7p5k_set3E60' 'ep2d_diff_2p5k_set4E60' };
            for ii=1:4 % 4 sets of DWIs in this project!
                obj.Params.DCM2NII.in(ii).prefix = dwi_setnames(ii); 
                [ ~ , obj.Params.DCM2NII.in(ii).nvols ] = system([ 'cat ' ...
                    obj.Params.DCM2NII.scanlog ' | grep ' dwi_setnames{ii} ...
                    ' | tail -1 | awk ''{ print $7 }'' ' ]);
                obj.Params.DCM2NII.in(ii).nvols=str2num(obj.Params.DCM2NII.in(ii).nvols);
                [ ~ , obj.Params.DCM2NII.in(ii).first_dcmfiles ]= system([ 'cat ' ...
                    obj.Params.DCM2NII.scanlog ' | grep ' dwi_setnames{ii} ...
                    ' | tail -1 | awk ''{ print $8 }'' ' ]);
                            
                obj.Params.DCM2NII.out(ii).location = [ obj.root 'Orig' filesep ];
                obj.Params.DCM2NII.out(ii).filename = [  cell2char(dwi_setnames(ii)) '.nii.gz' ];
            end
            
            %%%%%%%%%%%%
            %For proc_dropvols
            obj.Params.DropVols.in.dropVols=[ 1 ];
            obj.Params.DropVols.in.prefix='dv_';
            obj.Params.DropVols.in.movefiles=['..' filesep '01_DropVols' filesep ];
            obj.Params.DropVols.in.fn=obj.rawfiles;
            obj.Params.DropVols.out.fn=dir_wfp([obj.root, '01_DropVols', filesep, '*.nii.gz']);
            
            %%%%%%%%%%%%            
            %For gradient non-linearity correction
            obj.Params.GradNonlinCorrect.in.movefiles = '../02_GradCorrect/';
            obj.Params.GradNonlinCorrect.in.prefix = 'gnc_';
            obj.Params.GradNonlinCorrect.in.gradfile = obj.gradfile;
            
            obj.Params.GradNonlinCorrect.in.fslroi = [ 0 1 ];
            obj.Params.GradNonlinCorrect.in.target = obj.Params.DropVols.out.fn;
      
            %For proc_runFS
            obj.fsdir='/autofs/eris/bang/ADRC/FreeSurfer6.0/';
            obj.fsubj=obj.sessionname;
                
        end
              
        function obj = CommonProc(obj)
            
            %DICOM TO NIFTII
            obj.proc_dcm2nii;
            
            %DROPPING THE FIRST B0 VOLUME
            obj.proc_drop_vols(obj.rawfiles);
            
            %APPLYING NONLINEARITY CORRECTION (
            obj.proc_gradient_nonlin_correct
        end
        
       
        
        
        
        function resave(obj)
            save([obj.objectHome filesep obj.sessionname '.mat'],'obj');
        end
    end
end