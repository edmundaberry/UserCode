import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo2")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )

process.load("Configuration.StandardSequences.GeometryPilot2_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'CRAFT_V3P::All'

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

        '/store/data/Commissioning08/Calo/RECO/v1/000/068/288/009FCAA1-44A8-DD11-AF40-001617C3B6DC.root',
        '/store/data/Commissioning08/Calo/RECO/v1/000/068/288/00A4ED8C-3BA8-DD11-8E5F-000423D98834.root',
        '/store/data/Commissioning08/Calo/RECO/v1/000/068/288/00BBCE07-4BA8-DD11-AF54-001617E30D12.root',
        '/store/data/Commissioning08/Calo/RECO/v1/000/068/288/24860E46-56A8-DD11-B9EF-000423D98920.root',
        '/store/data/Commissioning08/Calo/RECO/v1/000/068/288/FC1FBA43-4FA8-DD11-9736-000423D98634.root',
        '/store/data/Commissioning08/Calo/RECO/v1/000/068/288/FEF4D585-42A8-DD11-96E1-000423D98E30.root'

        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_1.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_2.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_3.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_4.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_5.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_6.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_7.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_8.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_9.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_10.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_11.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_12.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_13.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_14.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_15.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_16.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_17.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_18.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_19.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_20.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_21.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_22.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_23.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_24.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_25.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_26.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_27.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_28.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_29.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_30.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_31.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_32.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_33.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_34.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_35.root',
        ##'dcache:/pnfs/cms/WAX/2/eberry/L1Skims/L1SkimRaw_ppex_36.root'
        ##'file:///uscms/home/eberry/scratch/L1SkimRaw_ppex_36.root'
    )
)

## Defines HLT reconstruction paths
process.load("HLTrigger.HLTanalyzers.HLTopen_cff")

process.demo = cms.EDAnalyzer('L1SkimAnalyzer')

process.analyze = cms.Path(process.demo)

process.schedule = cms.Schedule(process.DoHLTJets,process.analyze)
##process.schedule = cms.Schedule(process.analyze)
