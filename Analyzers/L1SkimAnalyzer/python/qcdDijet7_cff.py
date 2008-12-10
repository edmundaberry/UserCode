import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",skipEvents = cms.untracked.uint32(0),fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/0C89F526-6997-DD11-BD9D-0015C5E5B22E.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/12955449-D097-DD11-A5DC-0015C5E5B335.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/16724939-5697-DD11-BBD7-001D09677B79.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/1EFA61BF-4A97-DD11-876A-001D0967D6B1.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/26367013-6897-DD11-B26C-00E08133C424.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/38F9F456-CA97-DD11-BAB6-0015C5E5B335.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/3ABE9D8E-7397-DD11-8DD8-0019B9F4055B.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/4247329D-6297-DD11-8744-001D0967DAF3.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/4662D50F-6F97-DD11-A51D-0019B9F4055B.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/46B6C964-F997-DD11-9D28-001D09645F4A.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/526B5B90-7B97-DD11-94AE-001D0964612C.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/5819305C-7797-DD11-B756-00E081338BAA.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/58DB3166-5E97-DD11-A6E0-0015C5E673CC.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/6A0AD237-5697-DD11-8C8E-001D0967D9CC.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/6A1F021C-6997-DD11-B74E-0015C5E59E7F.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/6E36885D-7797-DD11-B980-0019B9F4055B.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/6EC8AF90-7B97-DD11-A11A-001D09645F4A.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/6ED1E5DE-7F97-DD11-B63E-0015C5EBA8FA.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/70D07F94-6697-DD11-9CC7-0015C5E5B9C5.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/847EA866-5E97-DD11-A3DD-0015C5E5B22E.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/88AF16DE-6197-DD11-A676-001D09645B19.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/8ACABF5C-7797-DD11-BF4C-0019B9F4055B.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/946E630A-5A97-DD11-B448-0015C5E5B9C5.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/960290BA-7397-DD11-AADA-001D09646049.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/983CC913-6997-DD11-B0CE-0015C5E5B9C5.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/9CC38DDC-DC97-DD11-976B-0015C5E5B335.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/9E91710B-5A97-DD11-A556-0015C5E5B288.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/A054D36C-5E97-DD11-8DF3-00093D1440BA.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/A21DA399-7397-DD11-929E-0015C5E5B9A7.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/ACC40159-7797-DD11-BDD7-001D09645F4A.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/AE0E4192-7B97-DD11-82DF-001D09645ABF.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/AEE21A95-6697-DD11-9277-0015C5E5B9C5.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/B21199D9-F297-DD11-9595-001D09645F4A.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/B8402E26-6997-DD11-98E4-0015C5E5B22E.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/BEF0888F-7B97-DD11-A520-001EC9AF95C7.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/C023A31D-6797-DD11-8E64-001D0968F760.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/C8ACD101-F097-DD11-BCA6-001D0964612C.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/D4E96CCC-DC97-DD11-9585-0015C5E5B288.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/D636D9D7-F297-DD11-8F7C-001D09645ABF.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/D6A00158-7797-DD11-9B83-001D09645A9D.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/DA807A82-C797-DD11-A33F-001D0967C8D3.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/DE43DE87-D897-DD11-972B-0019B9E4897D.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/E0CBA015-6897-DD11-A4E8-00E081344B92.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/E8E03890-D897-DD11-ACDF-0019B9E4878A.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/F0628264-5E97-DD11-A8A7-0015C5E5B9C5.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/F2870C59-7797-DD11-8427-001D0964612C.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/F2D99E19-DC97-DD11-AE8F-001D0964612C.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/FAF7B238-7497-DD11-A430-0019B9E487D7.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0000/FCC75F8F-7B97-DD11-B57F-001D09645A9D.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0001/0C4F95E7-5E99-DD11-AABF-0015C5E5B288.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0001/302B3841-6299-DD11-931C-0015C5E5B335.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0001/3EC5F85F-6799-DD11-ADFB-0015C5E5B288.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0001/4AA520EB-5E99-DD11-9A49-0015C5EC47A2.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0001/8A0823EB-5E99-DD11-8224-0015C5EC47A2.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0001/DA35B037-6299-DD11-8087-0030487BB7E8.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0001/E657957A-7599-DD11-A248-0015C5EC47A2.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0001/FC4D567B-7599-DD11-9E2C-0015C5E5B335.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/00D9EF1B-0294-DD11-A7F3-0015C5E5B22E.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/0AA770DC-E893-DD11-BB55-0015C5E5B22E.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/0EFE6D0B-0294-DD11-B23C-0015C5E5B22E.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/14188DD1-E893-DD11-9F36-0015C5E5B22E.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/16F7D20E-F893-DD11-8701-0030487C1174.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/221A0823-0294-DD11-B5A0-0015C5E59F84.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/241D010C-F893-DD11-92DF-0030487C117A.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/24E2DD0F-F893-DD11-B240-0030487BB550.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/2E1B95CF-E893-DD11-B08D-0015C5E5B22E.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/3838230B-F893-DD11-91D9-0030487BB542.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/3CB059A7-EA93-DD11-9E36-0015C5E673CC.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/484A630D-F893-DD11-AAE4-0030487C11A0.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/489EC8E0-E893-DD11-AB9C-0015C5E59F84.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/52FB00F1-E893-DD11-A793-0015C5E673CC.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/66DEF6EF-E893-DD11-AFDC-0015C5E59E7F.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/74B00DE0-E893-DD11-827F-0015C5E59F84.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/845D6138-0294-DD11-AE38-0030487C1154.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/8AABC426-0294-DD11-82A3-0015C5E673CC.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/9006A210-F893-DD11-8389-0030487BB550.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/907A0C23-0294-DD11-97A1-0015C5E59F84.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/C076520B-0294-DD11-A9BB-0015C5E5B22E.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/C6D6072C-0294-DD11-8F6B-0030487BB7E4.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/DA176861-F993-DD11-A5BF-00215A452926.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/DC8BEA28-0294-DD11-8F6B-0015C5E673CC.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/E417960C-F893-DD11-BBD0-0030487C1380.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/EC36FD0A-0294-DD11-9C11-0015C5E59E7F.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/F00A1F6C-0194-DD11-B0B7-0015C5E5B335.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/F09C030E-F893-DD11-8291-0030487C11A0.root',
       '/store/mc/Summer08/QCDDiJetPt120to170/GEN-SIM-RAW/IDEAL_V9_v1/0010/FA4AEC3A-0294-DD11-857C-0030487C1242.root' ] );


secFiles.extend( [
               ] )

