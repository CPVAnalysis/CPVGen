import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *

generator = cms.EDFilter(
    "Pythia8GeneratorFilter",
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    comEnergy = cms.double(13000.0), #TODO make sure to update for Run 3
    ExternalDecays = cms.PSet(
        EvtGen130 = cms.untracked.PSet(
            decay_table = cms.string('GeneratorInterface/EvtGenInterface/data/DECAY_2014_NOLONGLIFE.DEC'),
            particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt_2014.pdl'),
            #particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt_2020.pdl'), # available with 10_6_X
            user_decay_embedded= cms.vstring(
                # PVV_CPLH beta eta |G1+| argG1+ |G0+| argG0+ |G1-| argG1- 
                # beta: CKM angle 
                'Define Hp 0.49',
                'Define Hz 0.775',
                'Define Hm 0.4',
                'Define pHp 2.50',
                'Define pHz 0.0',
                'Define pHm -0.17',
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
            'PTFilter:scaleToFilter = 1.0' # in BHNL: is 5.0
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

# ask for the presence of a triggering muon in the event
#TODO enforce that this muon comes from the other B meson?
TriggerMuonFilter = cms.EDFilter("PythiaFilterMultiMother", 
    MaxEta = cms.untracked.double(1.55),
    MinEta = cms.untracked.double(-1.55),
    MinPt = cms.untracked.double(6.8), 
    ParticleID = cms.untracked.int32(13),
    MotherIDs = cms.untracked.vint32(511, 521, 531, 5122, 5112, 5212, 5222, 5132, 5232, 5332), # B mesons and B baryons
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

PhiDecayFilter = cms.EDFilter( #FIXME does not seem to be applied on both phi candidates
    "PythiaDauVFilter",
    verbose = cms.untracked.int32(0),
    NumberDaughters = cms.untracked.int32(2),
    MotherID = cms.untracked.int32(531),
    ParticleID = cms.untracked.int32(333),
    DaughterIDs = cms.untracked.vint32(321, -321), # kaons
    MinPt = cms.untracked.vdouble(0.5, 0.5), 
    MinEta = cms.untracked.vdouble(-2.5, -2.5),
    MaxEta = cms.untracked.vdouble(2.5, 2.5)
)

ProductionFilterSequence = cms.Sequence(generator*BFilter*BsDecayFilter*PhiDecayFilter*TriggerMuonFilter)
