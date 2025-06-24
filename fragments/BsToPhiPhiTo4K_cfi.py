import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *

generator = cms.EDFilter(
    "Pythia8GeneratorFilter",
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    comEnergy = cms.double(13000.0),
    ExternalDecays = cms.PSet(
        EvtGen130 = cms.untracked.PSet(
            decay_table = cms.string('GeneratorInterface/EvtGenInterface/data/DECAY_2014_NOLONGLIFE.DEC'),
            particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt_2020.pdl'),
            user_decay_embedded= cms.vstring(
                # The following relationshop must hold: 
                # HpHp + HzHz + HmHm = 1
                #
                # Latest values taken from PDG2024
                #   Amplitudes
                #     GammaL / Gamma = |A0|^2 = 0.379
                #     Gamma_perp / Gamma = |A_perp|^2 = 0.310
                #     Gamma_par = |A_par|^2 = 1 - GammaL - Gamma_perp = 0.311
                #     
                #     Definition: the amplitudes correspond to the transversity amplitudes rather than the helicity 
                #                 i.e. Hp = A_par, Hz = A0, Hm = A_perp
                #
                #                 Hp = sqrt(0.311) = 0.558
                #                 Hz = sqrt(0.379) = 0.616
                #                 Hm = sqrt(0.310) = 0.557
                #
                #   Phases
                #     phi_par = delta_par - delta0 = 2.469  
                #     phi_perp = delta_perp - delta0 = 2.75
                #
                #     Definition: pHp = delta_par, pHz = delta0, pHm = delta_perp
                #                 By convention, delta0 = 0
                #                 pHp = 2.469
                #                 pHz = 0.0
                #                 pHp = 2.75
                'Define Hp 0.558',
                'Define Hz 0.616',
                'Define Hm 0.557',
                'Define pHp 2.469',
                'Define pHz 0.0',
                'Define pHm 2.75',
                '',
                'Alias      MyB_s0   B_s0',
                'Alias      Myanti-B_s0   anti-B_s0',
                'ChargeConj Myanti-B_s0   MyB_s0',
                'Alias      MyPhi phi',
                'ChargeConj MyPhi MyPhi',
                '',
                'Decay MyB_s0',
                '1.000       MyPhi      MyPhi    PVV_CPLH 0.02 1 Hp pHp Hz pHz Hm pHm;', # info on decay model: https://evtgen.hepforge.org/doc/models.html
                'Enddecay',
                'Decay Myanti-B_s0',
                '1.000       MyPhi      MyPhi    PVV_CPLH 0.02 1 Hp pHp Hz pHz Hm pHm;',
                'Enddecay',
                'Decay MyPhi',
                 '1.000      K+         K-       VSS;', # decay of a vector particle into two scalar particles
                'Enddecay',
                'End'
            ), 
            list_forced_decays = cms.vstring(
                'MyB_s0', 
                'Myanti-B_s0', 
                'MyPhi'
            ),
            operates_on_particles = cms.vint32(531, -531, 333), 
            convertPythiaCodes = cms.untracked.bool(False),
        ),
        parameterSets = cms.vstring('EvtGen130')
    ),
    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        processParameters = cms.vstring(
            "SoftQCD:nonDiffractive = on",
            'PTFilter:filter = on',
            'PTFilter:quarkToFilter = 5',
            'PTFilter:scaleToFilter = 5.0'
        ),
        parameterSets = cms.vstring(
            'pythia8CommonSettings',
            'pythia8CP5Settings',
            'processParameters',
        )
    )
)

BFilter = cms.EDFilter(
    "PythiaFilter",
    MaxEta = cms.untracked.double(10.),
    MinEta = cms.untracked.double(-10.),
    ParticleID = cms.untracked.int32(531)
)

BsDecayFilter = cms.EDFilter(
    "PythiaDauVFilter",
    verbose         = cms.untracked.int32(1),
    NumberDaughters = cms.untracked.int32(2),
    ParticleID      = cms.untracked.int32(531),
    DaughterIDs     = cms.untracked.vint32(333, 333),      # phi
    MinPt           = cms.untracked.vdouble(0.0, 0.0),     # no filter for now
    MinEta          = cms.untracked.vdouble(-10.0, -10.0), # no filter for now
    MaxEta          = cms.untracked.vdouble(10.0,  10.0)   # no filter for now
)

PhiDecayFilter = cms.EDFilter(
    "MCMultiParticleFilter",
    ParticleID = cms.vint32(321),
    MotherID = cms.untracked.vint32(333),
    EtaMax = cms.vdouble(2.5),
    PtMin = cms.vdouble(0.3), # keep the threshold low
    Status = cms.vint32(0),
    NumRequired = cms.int32(4),
    AcceptMore = cms.bool(True)
)

TriggerMuonFilter = cms.EDFilter(
    "PythiaFilter",
    MaxEta = cms.untracked.double(1.55),
    MinEta = cms.untracked.double(-1.55),
    MinPt = cms.untracked.double(6.8),
    ParticleID = cms.untracked.int32(13)
)

ProductionFilterSequence = cms.Sequence(generator * BFilter * BsDecayFilter * PhiDecayFilter * TriggerMuonFilter)
