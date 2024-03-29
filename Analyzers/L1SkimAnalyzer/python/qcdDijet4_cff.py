import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",skipEvents = cms.untracked.uint32(0),fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/00D0ECB1-A197-DD11-8FC7-001EC9B0871C.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/069D46B2-A197-DD11-8103-001D0964474D.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/10AD8F88-0E97-DD11-A92B-0015C5E59F84.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/12860B7D-1897-DD11-A9BD-0015C5E5B288.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/148F2ECA-A197-DD11-8D45-001D0964474D.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1AEDB012-DC97-DD11-8685-0019B9F4055B.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/1CAC44A3-A197-DD11-B5F9-001D096460D5.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/2046F0B3-A197-DD11-B4FD-001D096460D5.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/221846B2-A197-DD11-95B6-001D0964474D.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/247330CA-A197-DD11-AED3-001D0964474D.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/248085BB-1297-DD11-A851-0015C5EC47A2.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/302ADB90-0E97-DD11-A003-00093D1440BA.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/38980C8A-0E97-DD11-8561-0015C5E5B961.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/38A14088-0A97-DD11-B786-0015C5E5B335.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/3AA38DAD-FD96-DD11-A7EF-0015C5E5B335.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/468974BA-1297-DD11-98F1-0015C5E5B288.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/4C07104F-3B97-DD11-A788-0019B9F4055B.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/52F65FA3-A197-DD11-A99A-001D0964474D.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/584E048A-0E97-DD11-924A-0015C5E5B961.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/627EB6B4-A197-DD11-8571-001D0964474D.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/6440F1EC-A197-DD11-A273-001D0964474D.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/64C67CB8-A197-DD11-8587-001EC9B0871C.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/6652988D-0A97-DD11-A89C-00215A45F882.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/6E4E63E5-A197-DD11-937B-001D096460D5.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/76FA748A-0A97-DD11-8660-0030487C2174.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/7A836EBA-1297-DD11-BE57-0015C5E5B335.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/7EF6D6B9-A197-DD11-9F9F-0019B9F4054C.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/88F1C4BA-1297-DD11-BA99-0015C5EC47A2.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/8A55C4D9-A197-DD11-9B82-001D0964474D.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/8ACF39BB-1297-DD11-9132-0015C5EC47A2.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/8C0F558A-0E97-DD11-8191-0015C5E673BD.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/94DEC0A4-4297-DD11-BA7B-001D0964612C.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/9A0B999A-0A97-DD11-8F1E-002264064196.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/A65916BB-1297-DD11-8746-0015C5E5B9C5.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/A6C1EBB1-A197-DD11-B6FF-001EC9B0871C.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/A6CADCBA-1297-DD11-9299-0015C5E5B335.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/AA4ECA98-0A97-DD11-9318-00215A45F8B6.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/AAF509B6-A197-DD11-9831-001D096460D5.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/AE1C54A1-0A97-DD11-BC71-0015C5EC47A2.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/AE1C79BA-1297-DD11-B4C9-0015C5E5B335.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/AE4EA6BA-1297-DD11-B040-0015C5E5B288.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/B03F40DF-A197-DD11-A3FC-001D0964474D.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/B465F8A7-A197-DD11-AB1C-001EC9B0871C.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/BA77CF8A-E496-DD11-8573-00215A45F882.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/BE0D6719-0697-DD11-9AB2-0015C5E5B961.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/C4673E89-0E97-DD11-B8C9-0015C5EC47A2.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/C81B4BBB-1297-DD11-861C-0015C5E5B288.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/CA8623A1-0A97-DD11-99D2-0015C5E5B335.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/CCBAB2BA-1297-DD11-9966-0015C5EC47A2.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/CE58E4F4-A197-DD11-B455-001EC9B0871C.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/D254774D-1A98-DD11-860D-001D0967DA85.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/D2C074BA-1297-DD11-B856-0015C5E5B335.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/D4F2FFE3-A197-DD11-BCD9-001D0964474D.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/D6BBD84D-BC97-DD11-9183-0030487BB7E4.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/DA60A4E9-A197-DD11-BCDC-001EC9B0871C.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/E6016F87-0A97-DD11-87D0-0015C5E5B335.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/E623E379-0497-DD11-B21B-0019B9F4055B.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/EAD4D5D1-1897-DD11-B2C3-0015C5EC47A2.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/EC8A4389-0E97-DD11-A6B0-0015C5E5B961.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/F4B3C2AE-A197-DD11-A881-001D096460D5.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/F4F139E0-A197-DD11-B77F-001EC9B0871C.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0000/FE8A708B-0E97-DD11-983D-0015C5E673BD.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0001/00494FB3-5598-DD11-8AB4-00E081334BC6.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0001/062472EA-5198-DD11-8EDF-001D09645ABF.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0001/2496F6FA-5E98-DD11-B739-0019B9F4056A.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0001/2C4347A1-FD98-DD11-AD96-001A6478AC14.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0001/2C82F7E2-5998-DD11-9172-001D09645ABF.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0001/34C0B8E8-5198-DD11-8F4B-0019B9F4055B.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0001/405569B6-5598-DD11-991D-001D09645F4A.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0001/4478BF77-AD98-DD11-9F57-001D09645A9D.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0001/48AEC156-BA98-DD11-89B9-001D09645B46.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0001/6CF5EB79-7598-DD11-91A3-00E081404000.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0001/74BFC14E-6298-DD11-8156-0019B9F4056A.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0001/7A958A0A-5298-DD11-A3BF-001D09645ABF.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0001/96A7C6FC-5E98-DD11-93E7-001D09645B46.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0001/9E4D0351-4D98-DD11-8D5B-001D09645B46.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0001/AA859832-3D98-DD11-BC3E-0015C5EB87DE.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0001/E4DFEF71-4998-DD11-B07E-001D09645F27.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0001/E8A5DE7D-B698-DD11-A743-001EC9AF95C7.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/0225D6D5-4792-DD11-AEFA-003048C17FC8.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/04F8EEE8-4692-DD11-8DD6-003048C17FC6.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/04FBC048-E992-DD11-888E-001D0967DAE4.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/0689A055-6392-DD11-81F9-003048C17FBC.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/06EF3E1C-1F93-DD11-A2EA-001D0967C130.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/08860171-4792-DD11-8AE2-003048C17FC6.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/0A69DD28-4792-DD11-8916-003048C17FC6.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/0AC326FD-2D93-DD11-AA25-001D096908D8.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/0C6A06D8-8F91-DD11-B6E4-00093D127D2D.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/0E99C63C-5C92-DD11-A3EA-003048C185DC.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/0EB23D03-F692-DD11-A9C6-0019B9E489B4.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/121C417C-A891-DD11-AD9E-0030487C1380.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/12D036E5-8B91-DD11-A93C-003048C26CB8.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/1657848D-9C91-DD11-8380-001E0B5A6378.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/1AD25F81-0193-DD11-BE02-0019B9E4A9BB.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/1E0355C4-7C91-DD11-8D10-003048C180D8.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/1EEE9611-5C92-DD11-90CD-003048C1872C.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/20860745-A191-DD11-B9CB-001E0B5A6378.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/228C6783-F492-DD11-AA7D-001D0967DE4A.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/2678BDE8-8B91-DD11-9D5E-003048C180D8.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/26C2F64C-4992-DD11-9FF7-0030487A9EA8.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/2E4339CF-5592-DD11-85AF-0030487C116E.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/2E6FAD46-4892-DD11-89EF-00304879C1EA.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/30369468-8191-DD11-87FE-003048C26CB6.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/309D9595-9191-DD11-8F3C-0015C5E5B288.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/32342705-0093-DD11-B040-0019B9E4A9BB.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/36AC47A7-9091-DD11-AA41-00093D127CD0.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/38DD7BD2-8F91-DD11-A5C8-00093D128C99.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/3AC70F56-6392-DD11-8324-003048C17FBC.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/40E8B98C-9391-DD11-BDC3-00093D128C99.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/4263E3A5-EB92-DD11-902B-001D096B0C83.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/443DE5CF-E892-DD11-951A-0019B9E49554.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/44436676-A891-DD11-8B7E-0030487BB542.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/4614EC57-7D91-DD11-A226-003048C18014.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/46F547CB-4992-DD11-AF72-0030487A9C52.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/48B9E6E8-8B91-DD11-B3D0-003048C180D8.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/4A90822A-5C92-DD11-8B87-003048C26CB8.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/4C6D0D54-8D91-DD11-B986-00304879F25A.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/4CDB012D-C591-DD11-B4F7-0015C5E5B335.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/52490E53-B091-DD11-8C64-00215A490902.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/54015711-8B91-DD11-B19B-003048C26CB8.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/58306154-8D91-DD11-9EE7-00304879F25A.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/58BC2952-4E92-DD11-9535-00093D13BB43.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/5EA343CA-0A93-DD11-9C9B-0019B9E4ABA6.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/5EFD6E7E-5E92-DD11-8A0F-0015C5E5B335.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/623902EE-DF91-DD11-B423-00215A4509DE.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/625F5E1D-F692-DD11-952D-001D0967D567.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/64011393-0493-DD11-B2B7-001D0967E002.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/6A59777B-0493-DD11-946E-001D0968F0B2.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/6C0A53F2-E191-DD11-B2F1-0030487C116E.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/6C0A7896-9191-DD11-9CA0-0015C5EC47A2.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/6CCFF4CB-7E91-DD11-B030-003048C180D8.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/74D96ABC-EB92-DD11-A604-001125C4617A.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/76D7FF62-8991-DD11-8BBB-003048C185DC.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/7880041D-5192-DD11-AB3C-001E0B5A5388.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/78F61073-E692-DD11-B6C2-0019B9E585EE.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/803AAD09-9691-DD11-A88B-00093D128C99.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/82F41DC3-4592-DD11-A862-003048C17FBC.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/8666728C-7B91-DD11-A287-003048C18014.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/86EB5EFD-5B92-DD11-8FAB-003048C185DC.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/86FF9A96-9191-DD11-A63D-0015C5EC47A2.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/8C6608BE-9E91-DD11-9CF5-001E0B5A6378.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/8EAA1A43-F292-DD11-B6DD-001D0967D238.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/92C6AF64-1693-DD11-9D33-001D0967DEA4.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/945338F9-4C92-DD11-99FB-00093D1440BA.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/94AA6979-0C93-DD11-9968-0019B9E48B41.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/9609168F-B891-DD11-AEB6-0030487BB7E8.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/9A782480-1C93-DD11-81FC-001D0967C103.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/9C8391F7-0793-DD11-A18E-0019B9E4B150.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/9CBD9E20-A391-DD11-A0A5-001E0B5A6378.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/A0761549-5092-DD11-8CCE-00093D145785.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/A2E2598B-9391-DD11-B635-00093D127CD0.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/A603B5AE-9091-DD11-B107-00093D127C69.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/A61AFA54-BB91-DD11-9398-0015C5E5B335.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/A8F16578-A891-DD11-B9F1-0030487C216E.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/AA26357C-DA91-DD11-A8E7-00093D1455B5.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/AAE39C60-C291-DD11-B0AA-0015C5E5B335.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/AC658C96-9191-DD11-BAB7-0015C5EC47A2.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/AC7E220C-9591-DD11-AC1A-00093D128C99.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/AE0FAE96-9191-DD11-AE67-0015C5E5B335.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/AE2EFAAE-B792-DD11-93AC-00192165CCD8.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/B070C3F9-F592-DD11-9C70-001D0967D021.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/B27A8549-D491-DD11-B02B-0030487C116E.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/B2D8212D-B091-DD11-9E3A-0030487C2174.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/B6736C6A-D491-DD11-8526-0030487BB7E4.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/B6845F3C-E992-DD11-A43B-001D0967DC6F.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/B816C97B-FC92-DD11-8E04-0019B9E7133E.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/B87DE5A6-9091-DD11-A3D0-00093D128C99.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/B8A5833C-5C92-DD11-BBDA-003048C185DC.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/BC1D420D-0093-DD11-AF3C-0019B9E48D5F.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/BC5B77A8-0093-DD11-8B74-001D0967B82E.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/BECE7088-4E92-DD11-AB64-00093D144DD2.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/CA89FDAF-EB92-DD11-94DD-0019B9E48CDB.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/CAF32AD1-F192-DD11-A4FB-001125C44CEA.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/CE6125E1-9491-DD11-A081-0030487A9EA4.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/D290FDA7-9091-DD11-96CF-00093D127C69.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/D441D37C-0493-DD11-BA3C-0019B9E4AC05.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/D873A494-8791-DD11-938D-003048C26CB8.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/D87D11EF-8891-DD11-B206-003048C185DC.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/D88185FE-0493-DD11-860D-0019B9E48D55.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/DA8B9242-4992-DD11-ABA9-0030487A9EA8.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/DA9C7A04-4892-DD11-891B-00304879C1EA.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/DEA99390-FF92-DD11-8AA2-0019B9E48DA0.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/E670A07C-F692-DD11-BF25-001D0967DE45.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/E6A3C0C4-F192-DD11-B37D-0019B9E48FC0.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/E8CF23B8-DC91-DD11-8D8A-0030487DE130.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/EC0BE644-F292-DD11-9591-001D0967D035.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/EEF95E18-A391-DD11-AB30-003048C26CB6.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/F001846C-7A91-DD11-AC23-003048C18014.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/F22F0B9F-FD92-DD11-A27B-0019B9E71483.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/F2788008-9691-DD11-849D-00093D127CD0.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/F6617E04-F492-DD11-8BC5-0019B9E4B05B.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/F66C078D-F592-DD11-BED0-0019B9E4B0EC.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/FAC952D0-0A93-DD11-BD27-0019B9E49111.root',
       '/store/mc/Summer08/QCDDiJetPt30to50/GEN-SIM-RAW/IDEAL_V9_v1/0009/FC488743-D491-DD11-A7C8-0030487BEF7A.root' ] );


secFiles.extend( [
               ] )

