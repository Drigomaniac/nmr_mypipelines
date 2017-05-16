classdef dwi_HAB < dwiMRI_Session
%%  classdef dwi_ADRC < dwiMRI_Session
%%  This class is a subclass of its parent class dwi_MRI_Session.m 
%%  (where it will inherent other methods). 
%%  Created by:
%%              Aaron Schultz aschultz@martinos.org
%%              Rodrigo Perea rpereacamargo@mgh.harvard.edu
%%
%%
    
    properties
       %root directoy where raw data lives:
       root = '/cluster/sperling/HAB/Project1/DWIs_30b700/Sessions/';
       dcm_location = '/cluster/sperling/HAB/Project1/DICOM_ARCHIVE/All_DICOMS/';
       session_location='/cluster/sperling/HAB/Project1/Sessions/'
    end
    methods
        function obj = dwi_HAB(sessionname,opt)
            %%%  If opt is passed, then the root Sessions folder will be 
            %%%  replaced with this argument.
            if nargin>1
                obj.root = opt;
            end
            
            obj.sessionname = sessionname;
            obj.root = [obj.root sessionname '/DWIs/'];
            obj.dcm_location= [ obj.dcm_location sessionname filesep ];
            obj.session_location= [ obj.session_location sessionname filesep ];
            %If the folder /DWIs/ does not exist, then create it!             
            if exist(obj.root,'dir')==0
                try
                    system(['mkdir -p' obj.root ]); 
                catch
                    disp([ 'Trying to create /DWIs/ in:' obj.root ... 
                        ' Maybe some permission issues?'])
                end
            end
            obj.objectHome = obj.root;
            
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
            %Global parameters:
            obj.vox= [2 2 2 ];
            obj.dosave=true;
            obj.rawfiles = dir_wfp([obj.root 'Orig/*.nii.gz' ] );
            if isempty(obj.rawfiles)
                obj.rawfiles = dir_wfp([obj.root 'Orig/*.nii' ] );
            end
            
            
            %For proc_DCM2NII:
            obj.Params.DCM2NII.specific_vols=35;
            obj.Params.DCM2NII.scanlog = [ obj.session_location  filesep 'LogFiles' ...
                filesep 'scan.log' ] ;
            if ~exist(obj.Params.DCM2NII.scanlog,'file')
                error(['No scanlog found in:' obj.Params.DCM2NII.scanlog '. Exiting...']);
            end
            dwi_setnames={ 'DIFFUSION_HighRes_30' };
            for ii=1:1 % Only 1 set on this sequence acquisiton (for HAB1)!
                obj.Params.DCM2NII.in(ii).prefix = dwi_setnames(ii); 
                try
                [ ~ , obj.Params.DCM2NII.in(ii).nvols ] = system([ 'cat ' ...
                    obj.Params.DCM2NII.scanlog ... 
                    ' | grep ' dwi_setnames{ii} ' | grep " 35 " | tail -1 | awk ''{ print $7 }'' ' ]);
                catch
                    disp('There is no DCMfile associated that contains 35 volumes .');
                    disp(['The specific session is:' obj.sessionname ]); 
                end
                obj.Params.DCM2NII.in(ii).nvols=str2num(obj.Params.DCM2NII.in(ii).nvols);
                [ ~ , obj.Params.DCM2NII.in(ii).first_dcmfiles ] = system([ 'cat ' ...
                    obj.Params.DCM2NII.scanlog ... 
                    ' | grep ' dwi_setnames{ii} ' | grep " 35 " | tail -1 | awk ''{ print $8 }'' ' ]);
                            
                obj.Params.DCM2NII.out(ii).location = [ obj.root 'Orig' filesep ];
                obj.Params.DCM2NII.out(ii).filename = [ cell2char(dwi_setnames(ii)) '.nii.gz' ];
            end
            
            %For proc_drop_vols
            %%%%THIS IS NOT IMPLEMENTED AS PART OF THE DIFFUSION PIPELINE
            %%%%FOR HAB1 project! (rdp20 dated 04/11/2017)
            
            
            %For proc_<XYZ?>
                      obj.fsdir='/autofs/eris/bang/ADRC/FreeSurfer6.0/';
            obj.fsubj=obj.sessionname;
                
        end
              
        function obj = CommonProc(obj)
            %%%
            obj.proc_t1_spm;
            %%%
            obj.Params.DropVols.in.dropVols = 1:10;
            obj.Params.DropVols.in.movefiles = '../01_DropVols/';
            obj.proc_drop_vols(obj.rawfiles);
            %%%
            obj.Params.Realign.in.movefiles = '../02_Realign/';
            obj.proc_realign(obj.Params.DropVols.out.fn);
            %%%
            obj.Params.Reslice.in.movefiles = '../03_Resliced/';
            obj.proc_reslice(obj.Params.DropVols.out.fn);
            %%%
            obj.Params.GradNonlinCorrect.in.movefiles = '../04_GradCorrect/';
            obj.Params.GradNonlinCorrect.in.prefix = 'gnc_';
            obj.Params.GradNonlinCorrect.in.gradfile = '/autofs/space/kant_004/users/ConnectomeScanner/Scripts/adrc_diff_prep/bash/gradient_nonlin_unwarp/gradient_coil_files/coeff_AS302.grad';
            obj.Params.GradNonlinCorrect.in.fn = obj.Params.Reslice.out.fn;
            obj.Params.GradNonlinCorrect.in.target = obj.Params.Reslice.out.meanimage;
            
            obj.Params.GradNonlinCorrect.out.warpfile = [];
            obj.Params.GradNonlinCorrect.out.meannii = [];
            obj.Params.GradNonlinCorrect.out.fn = [];
            obj.proc_gradient_nonlin_correct;
            %%%
            obj.Params.Implicit_Unwarp.in.movefiles = '../05_ImpUnwarp/';
            % obj.proc_implict_unwarping(obj.Params.GradNonlinCorrect.out.meanimage,  '/autofs/space/kant_004/users/ConnectomeScanner/Sessions/150430_8CS00315/restingState/05_ImpUnwarp/mmean.nii');
            obj.proc_implict_unwarping(obj.Params.GradNonlinCorrect.out.meanimage,  obj.Params.GradNonlinCorrect.out.fn);
            %%%
            obj.Params.Coreg.in.style = 'iuw';
            obj.Params.Coreg.in.movefiles = '../05_ImpUnwarp/';
            obj.proc_get_fs_labels;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.Params.ApplyNormNew.in.movefiles = '../06_Normed/';
            obj.Params.ApplyNormNew.in.regfile = obj.Params.spmT1_Proc.out.regfile;
            obj.Params.ApplyNormNew.in.prefix = 'nn2_';
            
            obj.Params.ApplyNormNew.in.fn = obj.Params.Implicit_Unwarp.out.newmean;
            obj.proc_applynorm_new;
            obj.Params.ApplyNormNew.out.normmean = obj.Params.ApplyNormNew.out.fn{1};
            
            obj.Params.ApplyNormNew.in.fn = obj.Params.Implicit_Unwarp.out.fn;
            obj.proc_applynorm_new;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.Params.Smooth.in.kernel = [6 6 6];
            obj.Params.Smooth.in.prefix = 'ss6_';
            obj.Params.Smooth.in.movefiles = '../07_Smoothed/';
            obj.proc_smooth(obj.Params.ApplyNormNew.out.fn);
            
%             obj.genTBRmaps(obj.Params.Smooth.out.fn);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             obj.Params.Smooth.in.kernel = [3 3 3];
%             obj.Params.Smooth.in.prefix = 'ss3_';
%             obj.Params.Smooth.in.movefiles = '../08_Smoothed/';
%             obj.proc_smooth(obj.Params.Implicit_Unwarp.out.fn);
            
            
%             obj.Params.ApplyReverseNormNew.in.movefiles = [obj.root '/08_Smoothed/templates/'];
%             obj.Params.ApplyReverseNormNew.in.fn = dir_wfp('/autofs/space/schopenhauer_002/users/MATLAB_Scripts/Atlas/fMRI/SchultzMaps/StandardTemplates/*.nii');
%             obj.Params.ApplyReverseNormNew.in.targ = [obj.Params.Smooth.out.fn{1} ',1'];
%             obj.Params.ApplyReverseNormNew.in.regfile = obj.Params.spmT1_Proc.out.iregfile;
%             obj.proc_apply_reservsenorm_new;
            
%             obj.genTBRmaps(obj.Params.Smooth.out.fn);

%             TBR(obj.Params.Smooth.out.fn,obj.Params.ApplyReverseNormNew.out.fn,[],[],'_Standard',[obj.root '/08_Smoothed/'],0);
%             TBR(obj.Params.Implicit_Unwarp.out.fn,obj.Params.ApplyReverseNormNew.out.fn,[],[],'_Standard',[obj.root '/05_ImpUnwarp/'],0);
    end   %COMMENTED BECAUSE ITS FMRI
        
        function obj = proc_gradient_nonlin_correct(obj)
            wasRun = false;
            target = obj.Params.GradNonlinCorrect.in.target;
            [m h] = openIMG(target); if h.pinfo(1)~=1; h.pinfo(1)=1; spm_write_vol(h,m); end
            [a b c] = fileparts(target);
            outpath = obj.getPath(a,obj.Params.GradNonlinCorrect.in.movefiles);

            
            % addpath(genpath('/autofs/space/kant_002/users/rperea/DrigoScripts/adrc_diff/adrc_diff_prep/'));
            % mris_gradient_nonlin__unwarp_volume__batchmode_HCPS_v3(target, [outpath 'gc_mean.nii'], 'coeff_AS302.grad');
            
            infile = target;
            outfile = [outpath 'gnc_' b c];
            gradfile = obj.Params.GradNonlinCorrect.in.gradfile;
            
            %%% Compute the grdient nonlinearity correction
            if exist([outpath b '_deform_grad_rel.nii'],'file')==0
                cmd=['sh /autofs/space/kant_004/users/ConnectomeScanner/Scripts/adrc_diff_prep/run_mris_gradient_nonlin__unwarp_volume__batchmode_ADRC_v3.sh ' ...
                    '/usr/pubsw/common/matlab/8.5 ' ...
                    infile ' ' outfile ' ' gradfile ' '];
                system(cmd);
                wasRun = true;
            end
            obj.Params.GradNonlinCorrect.out.warpfile = [outpath b '_deform_grad_rel.nii'];
            
            %%% Apply the correction to the mean image.
            if exist(outfile,'file')==0
                cmd = ['applywarp -i ' infile ' -r ' infile ' -o ' outfile ' -w ' obj.Params.GradNonlinCorrect.out.warpfile ' --interp=spline'];
                runFS(cmd,pwd,3);
                system(['gunzip ' outpath '*.gz']);
                wasRun = true;
            end
            obj.Params.GradNonlinCorrect.out.meanimage = outfile;
            
            %%% Apply correction to full dataset
            fn = obj.Params.Reslice.out.fn;
            for ii = 1:numel(fn);
                infile = fn{ii};
                [a b c] = fileparts(infile);
                outpath = obj.getPath(a,obj.Params.GradNonlinCorrect.in.movefiles);
                outfile = [outpath 'gnc_' b c];
                if exist(outfile,'file')==0
                    cmd = ['applywarp -i ' infile ' -r ' infile ' -o ' outfile ' -w ' obj.Params.GradNonlinCorrect.out.warpfile ' --interp=spline'];
                    runFS(cmd,pwd,3);
                    system(['gunzip ' outpath '*.gz']);
                    wasRun = true;
                end
                obj.Params.GradNonlinCorrect.out.fn{ii,1} = outfile;
            end
            
            obj.UpdateHist(obj.Params.GradNonlinCorrect,'proc_gradient_nonlin_correct',obj.Params.GradNonlinCorrect.out.warpfile,wasRun);
        end
        
        function resave(obj)
            save([obj.objectHome filesep obj.sessionname '.mat'],'obj');
        end
        
      
    end
end