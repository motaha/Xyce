#!/bin/sh

admsXml -e ../../xyceVersion.xml -e ../../xyceBasicTemplates.xml -e ../../xyceImplementationFile.xml psp103.va
admsXml -e ../../xyceVersion.xml -e ../../xyceBasicTemplates.xml -e ../../xyceHeaderFile.xml psp103.va
emacs N_DEV_ADMSPSP103VA.C --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
emacs N_DEV_ADMSPSP103VA.h --batch --eval="(require 'cc-mode)" --eval="(c-set-offset 'substatement-open 0)" --eval="(c-set-offset 'arglist-intro 3)" --eval="(c-set-offset 'innamespace -2)" --eval="(setq-default indent-tabs-mode nil)" --eval='(indent-region (point-min) (point-max) nil)' -f save-buffer
patch < make_PSP_usable.diff

# cp N_DEV_ADMSPSP103VA.h ../../../../src/DeviceModelPKG/ADMS/include/
# cp N_DEV_ADMSPSP103VA.C ../../../../src/DeviceModelPKG/ADMS/src/
