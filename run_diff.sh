echo PlGFfree
python bifrost_diff_tool.py calibrator-assign /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/PlGFfree/pilot_run1/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/PlGFfree/pilot_run1/output/results/
python bifrost_diff_tool.py calibrator-assign /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/PlGFfree/pilot_run1_run2/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/PlGFfree/pilot_run1_run2/output/results/
#python bifrost_diff_tool.py calibrator-assign /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/PlGFfree/pilot_run1_run2_final/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/PlGFfree/pilot_run1_run2_final/output/results/
python bifrost_diff_tool.py bridge-lot /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/PlGFfree/bridging/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/PlGFfree/bridging/output/results/

echo PlGFdiss
python bifrost_diff_tool.py calibrator-assign /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/PlGFdiss/pilot_run1_run2_final/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/PlGFdiss/pilot_run1_run2_final/output/results/
python bifrost_diff_tool.py bridge-lot /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/PlGFdiss/bridging/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/PlGFdiss/bridging/output/results/

echo sFlt
python bifrost_diff_tool.py calibrator-assign /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/sFlt/pilot_run1_final/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/sFlt/pilot_run1_final/output/results/
python bifrost_diff_tool.py bridge-lot /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/sFlt/bridging/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/sFlt/bridging/output/results/

echo KIM1
python bifrost_diff_tool.py calibrator-assign /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/KIM1/pilot_run1/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/KIM1/pilot_run1/output/results/
python bifrost_diff_tool.py calibrator-assign /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/KIM1/pilot_run1_run2/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/KIM1/pilot_run1_run2/output/results/
python bifrost_diff_tool.py calibrator-assign /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/KIM1/pilot_run1_run2_final/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/KIM1/pilot_run1_run2_final/output_deviation/results/
#python bifrost_diff_tool.py calibrator-assign /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/KIM1/pilot_run1_run2_final/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/KIM1/pilot_run1_run2_final/output/results/
python bifrost_diff_tool.py bridge-lot /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/KIM1/bridging/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/KIM1/bridging/output/results/

echo FGF21
python bifrost_diff_tool.py calibrator-assign /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/FGF21/pilot_run1/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/FGF21/pilot_run1/output/results/
python bifrost_diff_tool.py calibrator-assign /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/FGF21/pilot_run1_final/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/FGF21/pilot_run1_final/output/results/
python bifrost_diff_tool.py bridge-lot /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/FGF21/bridging/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/FGF21/bridging/output/results/

echo CD274
python bifrost_diff_tool.py calibrator-assign /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/CD274/pilot_run1/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/CD274/pilot_run1/output/results/
python bifrost_diff_tool.py calibrator-assign /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/CD274/pilot_run1_final/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/CD274/pilot_run1_final/output/results/
python bifrost_diff_tool.py bridge-lot /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/CD274/bridging/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/CD274/bridging/output/results/

echo ENG
python bifrost_diff_tool.py calibrator-assign /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/ENG/pilot_run1/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/ENG/pilot_run1/output/results/
#python bifrost_diff_tool.py calibrator-assign /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/ENG/pilot_run1_run2/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/ENG/pilot_run1_run2/output/results/
python bifrost_diff_tool.py calibrator-assign /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/ENG/pilot_run1_run2_final/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/ENG/pilot_run1_run2_final/output/results/
python bifrost_diff_tool.py bridge-lot /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/ENG/bridging/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/ENG/bridging/output/results/

echo DCN
python bifrost_diff_tool.py calibrator-assign /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/DCN/pilot_run1/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/DCN/pilot_run1/output/results/
python bifrost_diff_tool.py calibrator-assign /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/DCN/pilot_run1_run2/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/DCN/pilot_run1_run2/output/results/
python bifrost_diff_tool.py calibrator-assign /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/DCN/pilot_run1_run2_final/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/DCN/pilot_run1_run2_final/output/results/
python bifrost_diff_tool.py bridge-lot /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/DCN/bridging/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/DCN/bridging/output/results/


#python bifrost_diff_tool.py calibrator-assign /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/PlGFfree/pilot_run1_run2_final/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/PlGFfree/pilot_run1_run2_final/output/results/
# python bifrost_diff_tool.py calibrator-assign /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/ENG/pilot_run1_run2/expected/ /mnt/ruo_rw/rnd/SQA_Data/output/progenity_app_bifrost/ENG/pilot_run1_run2/output/results/
