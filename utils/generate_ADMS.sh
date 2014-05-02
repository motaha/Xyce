#/bin/sh

# Be sure to trun this from the project root!

(cd utils/ADMS/examples/fbh_hbt-2.1 && ./make_FBH_usable.sh)
(cd utils/ADMS/examples/vbic/vbic_code/vbic_3T_et_cf && ./make_typed_vbic_usable.sh)
(cd utils/ADMS/examples/bsimcmg_107.0.0/code && ./make_bsimcmg_usable.sh)
(cd utils/ADMS/examples/psp103 && ./make_PSP_usable.sh)
(cd src/DeviceModelPKG/Xyce_NonFree/Verilog/ekv301_02 && ./make_ekv_usable.sh)
(cd src/DeviceModelPKG/Xyce_NonFree/Verilog/ekv2.6 && ./make_ekv2.6_usable.sh)

diff -u -Bb utils/ADMS/examples/fbh_hbt-2.1/N_DEV_ADMSHBT_X.h src/DeviceModelPKG/ADMS/include/N_DEV_ADMSHBT_X.h
diff -u -Bb utils/ADMS/examples/fbh_hbt-2.1/N_DEV_ADMSHBT_X.C src/DeviceModelPKG/ADMS/src/N_DEV_ADMSHBT_X.C
diff -u -Bb utils/ADMS/examples/vbic/vbic_code/vbic_3T_et_cf/N_DEV_ADMSvbic.h src/DeviceModelPKG/ADMS/include/N_DEV_ADMSvbic.h
diff -u -Bb utils/ADMS/examples/vbic/vbic_code/vbic_3T_et_cf/N_DEV_ADMSvbic.C src/DeviceModelPKG/ADMS/src/N_DEV_ADMSvbic.C
diff -u -Bb utils/ADMS/examples/bsimcmg_107.0.0/code/N_DEV_ADMSbsimcmg.h src/DeviceModelPKG/ADMS/include/N_DEV_ADMSbsimcmg.h
diff -u -Bb utils/ADMS/examples/bsimcmg_107.0.0/code/N_DEV_ADMSbsimcmg.C src/DeviceModelPKG/ADMS/src/N_DEV_ADMSbsimcmg.C
diff -u -Bb utils/ADMS/examples/psp103/N_DEV_ADMSPSP103VA.h src/DeviceModelPKG/ADMS/include/N_DEV_ADMSPSP103VA.h
diff -u -Bb utils/ADMS/examples/psp103/N_DEV_ADMSPSP103VA.C src/DeviceModelPKG/ADMS/src/N_DEV_ADMSPSP103VA.C
diff -u -Bb src/DeviceModelPKG/Xyce_NonFree/Verilog/ekv301_02/N_DEV_ADMSekv3.h src/DeviceModelPKG/Xyce_NonFree/include/N_DEV_ADMSekv3.h
diff -u -Bb src/DeviceModelPKG/Xyce_NonFree/Verilog/ekv301_02/N_DEV_ADMSekv3.C src/DeviceModelPKG/Xyce_NonFree/src/N_DEV_ADMSekv3.C
diff -u -Bb  src/DeviceModelPKG/Xyce_NonFree/Verilog/ekv2.6/N_DEV_ADMSekv_va.h src/DeviceModelPKG/Xyce_NonFree/include/N_DEV_ADMSekv_va.h
diff -u -Bb  src/DeviceModelPKG/Xyce_NonFree/Verilog/ekv2.6/N_DEV_ADMSekv_va.C src/DeviceModelPKG/Xyce_NonFree/src/N_DEV_ADMSekv_va.C

exit

cp utils/ADMS/examples/fbh_hbt-2.1/N_DEV_ADMSHBT_X.h src/DeviceModelPKG/ADMS/include/
cp utils/ADMS/examples/fbh_hbt-2.1/N_DEV_ADMSHBT_X.C src/DeviceModelPKG/ADMS/src/
cp utils/ADMS/examples/vbic/vbic_code/vbic_3T_et_cf/N_DEV_ADMSvbic.h src/DeviceModelPKG/ADMS/include/
cp utils/ADMS/examples/vbic/vbic_code/vbic_3T_et_cf/N_DEV_ADMSvbic.C src/DeviceModelPKG/ADMS/src
cp utils/ADMS/examples/bsimcmg_107.0.0/code/N_DEV_ADMSbsimcmg.h src/DeviceModelPKG/ADMS/include/
cp utils/ADMS/examples/bsimcmg_107.0.0/code/N_DEV_ADMSbsimcmg.C src/DeviceModelPKG/ADMS/src/
cp utils/ADMS/examples/psp103/N_DEV_ADMSPSP103VA.h src/DeviceModelPKG/ADMS/include/
cp utils/ADMS/examples/psp103/N_DEV_ADMSPSP103VA.C src/DeviceModelPKG/ADMS/src/
cp src/DeviceModelPKG/Xyce_NonFree/Verilog/ekv301_02/N_DEV_ADMSekv3.h src/DeviceModelPKG/Xyce_NonFree/include/
cp src/DeviceModelPKG/Xyce_NonFree/Verilog/ekv301_02/N_DEV_ADMSekv3.C src/DeviceModelPKG/Xyce_NonFree/src/
