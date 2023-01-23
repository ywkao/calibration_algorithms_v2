import FWCore.ParameterSet.Config as cms
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import hgceeDigitizer, hgchefrontDigitizer, hgchebackDigitizer, hfnoseDigitizer

from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9
process = cms.Process('PROD',Phase2C11I13M9)

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D86Reco_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 100 )

path = "file:/eos/user/y/ykao/www/HGCAL_Geant4_project"
tag = "D86_R90To130_E100"

process.source = cms.Source("PoolSource",
        #fileNames = cms.untracked.vstring('file:/eos/user/y/ykao/www/HGCAL_Geant4_project/testbeam_positron_D86_R80To100_E100/step2.root')
        fileNames = cms.untracked.vstring( path + '/testbeam_positron_' + tag + '/step2.root')
        )

#proccess.source = cms.Source("EmptySource")

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1))

process.DummyRecHitAnalyzer = cms.EDAnalyzer('DummyRecHitAnalyzer',
                                        DataType = cms.string("beam"),
                                        CalibrationFlags = cms.vint32(1, 1, 0, 0, 0, 0, 0, 0, 0, 0)

                                        #---------------- Options of DataType -----------------#
                                        # option 1: beam
                                        # option 2: pedestal
                                        #
                                        #---------- Definitions of calibration flags ----------#
                                        # calibration_flags[0]: pedestal subtraction
                                        # calibration_flags[1]: cm subtraction
                                        # calibration_flags[2]: BX-1 correction
                                        # calibration_flags[3]: gain linearization
                                        # calibration_flags[4]: charge collection efficiency
                                        # calibration_flags[5]: MIP scale
                                        # calibration_flags[6]: EM scale
                                        # calibration_flags[7]: zero suppression
                                        # calibration_flags[8]: hit energy calibration
                                        # calibration_flags[9]: ToA conversion
                                        #------------------------------------------------------#
        )

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('./rootfiles/output_pure_test.root')
        )

process.p = cms.Path(process.DummyRecHitAnalyzer)
