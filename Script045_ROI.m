clear all
close all
path(path, '/mnt/data/scratch/igilab/jslin1/RadPath/Functions')
load_sql
ROI_dir = [script01_prefix '04_Masks/ROIs_'];
script04_results = [script04_prefix 'Results/' ];
if ~isdir(script04_results)
    mkdir(script04_results)
end 

for jj = 1: size(series_T2_register,1)  
    [ROI_ptno, ROI_Descrip, ROI_UID, ROI_ptno_Descrip_UID] = Generate_Label(series_T2_register(jj));
    ROIs = load_nii([ ROI_dir  ROI_ptno_Descrip_UID '.nii.Labels.nii' ]);
    idx1 = find(ROIs.img==1); % csf
    idx2 = find(ROIs.img==2); % gm
    idx3 = find(ROIs.img==3); % wm
    idx4 = find(ROIs.img==4); % fat
    
    lot=1;
    prefix=script01_prefix;
    [csf1, gm1, wm1, fat1] = ROI_sample(lot, prefix, series_T1,         ROI_ptno, idx1, idx2, idx3, idx4);
    [csf2, gm2, wm2, fat2] = ROI_sample(lot, prefix, series_T1post,     ROI_ptno, idx1, idx2, idx3, idx4);
    [csf3, gm3, wm3, fat3] = ROI_sample(0,   prefix, series_T2_register,ROI_ptno, idx1, idx2, idx3, idx4);
    [csf4, gm4, wm4, fat4] = ROI_sample(lot, prefix, series_T2star,     ROI_ptno, idx1, idx2, idx3, idx4);
    [csf5, gm5, wm5, fat5] = ROI_sample(lot, prefix, series_T2FLAIR,    ROI_ptno, idx1, idx2, idx3, idx4);
    [csf6, gm6, wm6, fat6] = ROI_sample(lot, prefix, series_SWAN,       ROI_ptno, idx1, idx2, idx3, idx4);

    lot=2;
    prefix=script02_prefix;
    [csf7, gm7, wm7, fat7] = ROI_sample(lot, prefix, series_ADC,        ROI_ptno, idx1, idx2, idx3, idx4);
    [csf8, gm8, wm8, fat8] = ROI_sample(lot, prefix, series_eADC,       ROI_ptno, idx1, idx2, idx3, idx4);
    [csf9, gm9, wm9, fat9] = ROI_sample(lot, prefix, series_FA,         ROI_ptno, idx1, idx2, idx3, idx4);
    [csf10, gm10, wm10, fat10] = ROI_sample(lot, prefix, series_AvgDC,  ROI_ptno, idx1, idx2, idx3, idx4);
    
    lot=3;
    prefix=script03_prefix;
    [csf11, gm11, wm11, fat11] = ROI_sample(lot, prefix, series_rBV,      ROI_ptno, idx1, idx2, idx3, idx4);
    [csf12, gm12, wm12, fat12] = ROI_sample(lot, prefix, series_rBV_corr, ROI_ptno, idx1, idx2, idx3, idx4);
    [csf13, gm13, wm13, fat13] = ROI_sample(lot, prefix, series_rBF,      ROI_ptno, idx1, idx2, idx3, idx4);
    [csf14, gm14, wm14, fat14] = ROI_sample(lot, prefix, series_MTT,      ROI_ptno, idx1, idx2, idx3, idx4);
    [csf15, gm15, wm15, fat15] = ROI_sample(lot, prefix, series_Delay,    ROI_ptno, idx1, idx2, idx3, idx4);
    [csf16, gm16, wm16, fat16] = ROI_sample(lot, prefix, series_K2,       ROI_ptno, idx1, idx2, idx3, idx4);
    
    lot=3;
    prefix=script03_prefix;
    [csf17, gm17, wm17, fat17] = ROI_sample(lot, prefix, series_Ktrans,  ROI_ptno, idx1, idx2, idx3, idx4);
    [csf18, gm18, wm18, fat18] = ROI_sample(lot, prefix, series_Kep,     ROI_ptno, idx1, idx2, idx3, idx4);
    [csf19, gm19, wm19, fat19] = ROI_sample(lot, prefix, series_Vp,      ROI_ptno, idx1, idx2, idx3, idx4);
    [csf20, gm20, wm20, fat20] = ROI_sample(lot, prefix, series_Ve,      ROI_ptno, idx1, idx2, idx3, idx4);
    [csf21, gm21, wm21, fat21] = ROI_sample(lot, prefix, series_Wash_in, ROI_ptno, idx1, idx2, idx3, idx4);
    [csf22, gm22, wm22, fat22] = ROI_sample(lot, prefix, series_Wash_out,ROI_ptno, idx1, idx2, idx3, idx4);
    [csf23, gm23, wm23, fat23] = ROI_sample(lot, prefix, series_TTP,     ROI_ptno, idx1, idx2, idx3, idx4);
    [csf24, gm24, wm24, fat24] = ROI_sample(lot, prefix, series_AUC,     ROI_ptno, idx1, idx2, idx3, idx4);
    [csf25, gm25, wm25, fat25] = ROI_sample(lot, prefix, series_Peak,    ROI_ptno, idx1, idx2, idx3, idx4);

    csf_total = [csf1;      csf2;       csf3;       csf4;       csf5;...
                  csf6;      csf7;       csf8;       csf9;       csf10;...
                  csf11;     csf12;      csf13;      csf14;      csf15;...
                  csf16;     csf17;      csf18;      csf19;      csf20;...
                  csf21;     csf22;      csf23;      csf24;      csf25];
    gm_total = [gm1;      gm2;       gm3;       gm4;       gm5;...
                 gm6;      gm7;       gm8;       gm9;       gm10;...
                 gm11;     gm12;      gm13;      gm14;      gm15;...
                 gm16;     gm17;      gm18;      gm19;      gm20;...
                 gm21;     gm22;      gm23;      gm24;      gm25];
    wm_total = [wm1;      wm2;       wm3;       wm4;       wm5;...
                 wm6;      wm7;       wm8;       wm9;       wm10;...
                 wm11;     wm12;      wm13;      wm14;      wm15;...
                 wm16;     wm17;      wm18;      wm19;      wm20;...
                 wm21;     wm22;      wm23;      wm24;      wm25];
    fat_total = [fat1;      fat2;       fat3;       fat4;      fat5;...
                  fat6;      fat7;       fat8;       fat9;       fat10;...
                  fat11;     fat12;      fat13;      fat14;      fat15;...
                  fat16;     fat17;      fat18;      fat19;      fat20;...
                  fat21;     fat22;      fat23;      fat24;      fat25];
end

save([script04_results 'ROIs_' img_ptno '.mat'], 'csf_total','gm_total','wm_total','fat_total')

