import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_1.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_2.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_3.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_4.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_5.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_6.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_7.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_8.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_9.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_10.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_11.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_12.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_13.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_14.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_15.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_16.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_17.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_18.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_19.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_20.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_21.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_22.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_23.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_24.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_25.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_26.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_27.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_28.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_29.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_30.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_31.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_32.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_33.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_34.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_35.root',
        'file:///mnt/cms/cms1/eberry/L1SkimRaw_ppex_36.root'
    )
)

process.load("CondCore.DBCommon.CondDBSetup_cfi")

process.PoolDBESSource = cms.ESSource("PoolDBESSource",
    process.CondDBSetup,
    toGet = cms.VPSet(
        cms.PSet(
            record = cms.string('L1JetEtScaleRcd'),
            tag = cms.string('L1JetEtScale_CRUZET_hlt'),
            type = cms.string('L1CaloEtScale')
            ),
        cms.PSet(
            record = cms.string('L1EmEtScaleRcd'),
            tag = cms.string('L1EmEtScale_CRUZET_hlt'),
            type = cms.string('L1CaloEtScale')
            )),
    connect = cms.string('frontier://cmsfrontier.cern.ch:8000/FrontierProd/CMS_COND_20X_L1T')
)

process.demo = cms.EDAnalyzer('L1SkimAnalyzer',
                              l1CaloEmCandsTag = cms.untracked.InputTag('gctDigis')
)



process.p = cms.Path(process.demo)
