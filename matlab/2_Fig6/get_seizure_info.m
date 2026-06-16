function seizures = get_seizure_info(rootDir)

% =========================================================
% 保持你原来的 caseInfo 风格
% 关键：这里的 pre_idx / ictal_idx / post_idx
% 是基于“原始参考窗长”得到的窗口索引
%
% 如果你原始分析就是：
%   - 窗长 = 3 秒
%   - 步长 = 1 秒
% 那就保持下面这两个参数不变
% =========================================================
refWindowSec = 3;
refStepSec   = 1;

caseInfo = struct([]);

% ===================== 例子 =====================
caseInfo(1).case_id   = 'seizure_01';
caseInfo(1).file_path = fullfile('data', 'seizure_01_Gamma.mat');
caseInfo(1).pre_idx   = 1:200;
caseInfo(1).ictal_idx = 201:400;
caseInfo(1).post_idx  = 401:598;

caseInfo(2).case_id   = 'seizure_02';
caseInfo(2).file_path = fullfile('data', 'seizure_02_Gamma.mat');
caseInfo(2).pre_idx   = 1:122;
caseInfo(2).ictal_idx = 123:233;
caseInfo(2).post_idx  = 234:338;

caseInfo(3).case_id   = 'seizure_03';
caseInfo(3).file_path = fullfile('data', 'seizure_03_Gamma.mat');
caseInfo(3).pre_idx   = 1:195;
caseInfo(3).ictal_idx = 196:330;
caseInfo(3).post_idx  = 331:477;

caseInfo(4).case_id   = 'seizure_04';
caseInfo(4).file_path = fullfile('data', 'seizure_04_Gamma.mat');
caseInfo(4).pre_idx   = 1:137;
caseInfo(4).ictal_idx = 138:297;
caseInfo(4).post_idx  = 298:398;

caseInfo(5).case_id   = 'seizure_05';
caseInfo(5).file_path = fullfile('data', 'seizure_05_Gamma.mat');
caseInfo(5).pre_idx   = 1:252;
caseInfo(5).ictal_idx = 253:413;
caseInfo(5).post_idx  = 414:498;

caseInfo(6).case_id   = 'seizure_06';
caseInfo(6).file_path = fullfile('data', 'seizure_06_Gamma.mat');
caseInfo(6).pre_idx   = 1:199;
caseInfo(6).ictal_idx = 200:285;
caseInfo(6).post_idx  = 286:398;

caseInfo(7).case_id   = 'seizure_07';
caseInfo(7).file_path = fullfile('data', 'seizure_07_Gamma.mat');
caseInfo(7).pre_idx   = 1:239;
caseInfo(7).ictal_idx = 240:331;
caseInfo(7).post_idx  = 332:498;

caseInfo(8).case_id   = 'seizure_08';
caseInfo(8).file_path = fullfile('data', 'seizure_08_Gamma.mat');
caseInfo(8).pre_idx   = 1:236;
caseInfo(8).ictal_idx = 237:313;
caseInfo(8).post_idx  = 314:498;

caseInfo(9).case_id   = 'seizure_09';
caseInfo(9).file_path = fullfile('data', 'seizure_09_Gamma.mat');
caseInfo(9).pre_idx   = 1:59;
caseInfo(9).ictal_idx = 60:112;
caseInfo(9).post_idx  = 113:385;

caseInfo(10).case_id   = 'seizure_10';
caseInfo(10).file_path = fullfile('data', 'seizure_10_Gamma.mat');
caseInfo(10).pre_idx   = 1:244;
caseInfo(10).ictal_idx = 245:470;
caseInfo(10).post_idx  = 471:598;

caseInfo(11).case_id   = 'seizure_11';
caseInfo(11).file_path = fullfile('data', 'seizure_11_Gamma.mat');
caseInfo(11).pre_idx   = 1:150;
caseInfo(11).ictal_idx = 151:271;
caseInfo(11).post_idx  = 272:404;

caseInfo(12).case_id   = 'seizure_12';
caseInfo(12).file_path = fullfile('data', 'seizure_12_Gamma.mat');
caseInfo(12).pre_idx   = 1:180;
caseInfo(12).ictal_idx = 181:261;
caseInfo(12).post_idx  = 262:498;

caseInfo(13).case_id   = 'seizure_13';
caseInfo(13).file_path = fullfile('data', 'seizure_13_Gamma.mat');
caseInfo(13).pre_idx   = 1:69;
caseInfo(13).ictal_idx = 70:130;
caseInfo(13).post_idx  = 131:198;

caseInfo(14).case_id   = 'seizure_14';
caseInfo(14).file_path = fullfile('data', 'seizure_14_Gamma.mat');
caseInfo(14).pre_idx   = 1:252;
caseInfo(14).ictal_idx = 253:346;
caseInfo(14).post_idx  = 347:448;

caseInfo(15).case_id   = 'seizure_15';
caseInfo(15).file_path = fullfile('data', 'seizure_15_Gamma.mat');
caseInfo(15).pre_idx   = 1:154;
caseInfo(15).ictal_idx = 155:230;
caseInfo(15).post_idx  = 231:324;

caseInfo(16).case_id   = 'seizure_16';
caseInfo(16).file_path = fullfile('data', 'seizure_16_Gamma.mat');
caseInfo(16).pre_idx   = 1:253;
caseInfo(16).ictal_idx = 254:331;
caseInfo(16).post_idx  = 332:489;

caseInfo(17).case_id   = 'seizure_17';
caseInfo(17).file_path = fullfile('data', 'seizure_17_Gamma.mat');
caseInfo(17).pre_idx   = 1:135;
caseInfo(17).ictal_idx = 136:270;
caseInfo(17).post_idx  = 271:309;

caseInfo(18).case_id   = 'seizure_18';
caseInfo(18).file_path = fullfile('data', 'seizure_18_Gamma.mat');
caseInfo(18).pre_idx   = 1:117;
caseInfo(18).ictal_idx = 118:300;
caseInfo(18).post_idx  = 301:598;

caseInfo(19).case_id   = 'seizure_19';
caseInfo(19).file_path = fullfile('data', 'seizure_19_Gamma.mat');
caseInfo(19).pre_idx   = 1:255;
caseInfo(19).ictal_idx = 256:482;
caseInfo(19).post_idx  = 483:598;

caseInfo(20).case_id   = 'seizure_20';
caseInfo(20).file_path = fullfile('data', 'seizure_20_Gamma.mat');
caseInfo(20).pre_idx   = 1:311;
caseInfo(20).ictal_idx = 312:540;
caseInfo(20).post_idx  = 541:598;

caseInfo(21).case_id   = 'seizure_21';
caseInfo(21).file_path = fullfile('data', 'seizure_21_Gamma.mat');
caseInfo(21).pre_idx   = 1:94;
caseInfo(21).ictal_idx = 95:170;
caseInfo(21).post_idx  = 171:318;

caseInfo(22).case_id   = 'seizure_22';
caseInfo(22).file_path = fullfile('data', 'seizure_22_Gamma.mat');
caseInfo(22).pre_idx   = 1:95;
caseInfo(22).ictal_idx = 96:195;
caseInfo(22).post_idx  = 196:350;

caseInfo(23).case_id   = 'seizure_23';
caseInfo(23).file_path = fullfile('data', 'seizure_23_Gamma.mat');
caseInfo(23).pre_idx   = 1:162;
caseInfo(23).ictal_idx = 163:228;
caseInfo(23).post_idx  = 229:310;

caseInfo(24).case_id   = 'seizure_24';
caseInfo(24).file_path = fullfile('data', 'seizure_24_Gamma.mat');
caseInfo(24).pre_idx   = 1:29;
caseInfo(24).ictal_idx = 30:95;
caseInfo(24).post_idx  = 96:148;

% ===================== 例子 =====================

if isempty(caseInfo)
    seizures = struct([]);
    return;
end

seizures = struct([]);
for k = 1:numel(caseInfo)

    seizures(k).id   = caseInfo(k).case_id;
    seizures(k).file = resolve_file_path(rootDir, caseInfo(k).file_path);

    seizures(k).pre_idx   = caseInfo(k).pre_idx;
    seizures(k).ictal_idx = caseInfo(k).ictal_idx;
    seizures(k).post_idx  = caseInfo(k).post_idx;

    % =====================================================
    % 用“原始参考窗长 + 原始窗口索引”换算出 ictal 的参考中心时间边界
    % 后面无论新窗长是 2/3/4/5 秒，都按新窗中心时间去比较
    % =====================================================
    seizures(k).ictalStartCenterSec = ((min(caseInfo(k).ictal_idx) - 1) * refStepSec) + refWindowSec/2;
    seizures(k).ictalEndCenterSec   = ((max(caseInfo(k).ictal_idx) - 1) * refStepSec) + refWindowSec/2;

    seizures(k).refWindowSec = refWindowSec;
    seizures(k).refStepSec   = refStepSec;
end

end


function filePath = resolve_file_path(rootDir, filePathIn)

if isfile(filePathIn)
    filePath = filePathIn;
    return;
end

filePath = fullfile(rootDir, filePathIn);

end