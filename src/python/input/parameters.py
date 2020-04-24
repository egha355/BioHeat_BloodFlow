class Problem_Params:
    def __init__(self):
        self.coupled = 3 # CoupledBioheatFlow
        self.coupled_set(self.coupled)
        self.tissueMeshNumber=2 # (1)
        self.numberOfElementNodes = 4
        self.number_of_element_nodes_set(self.tissueMeshNumber)

        # time parameters, bioheat
        self.timeIncrementBioheat = 100 # (2)
        self.startTimeBioheat = 0.0
        self.timeStepsBioheat = 61 # (3)
        self.diffusionOutputFrequency = 1 # (4)

        # time parameters, flow
        self.flowOutputFrequency = 100 # (5)
        self.startTimeFlow = 0.0
        self.timeStepsFlow = 1 # (6)
        self.timeIncrementFlow = 0.1 # (7)
       
        # input mesh files
        self.flowNodeFile = 'input/Flow/Node.csv' # (8)
        self.flowElementFile = 'input/Flow/Element.csv' # (9)

        # properties
        self.conductivity_blood   = 0.5
        self.rho_blood            = 1069.0
        self.c_blood              = 3650.0
        self.rho_muscle           = 1085.0   #   muscle density
        self.c_muscle             = 3768.0        # J/Kg.K   muscle specific heat
        self.rho_bone             = 1357.0   # kg/m3    bone density
        self.c_bone               = 1700.0        # J/Kg.K   bone specific heat
        self.rho_skin             = 1085.0   # kg/m3    skin density
        self.c_skin               = 3768.0        # J/Kg.K   skin specific heat

        self.k_muscle             = 0.42     # W/m.K muscle conductivity.
        self.k_bone               = 0.75     # W/m.K bone conductivity.
        self.k_skin               = 0.42     # W/m.K skin conductivity.

        self.h_conv             = 100000.0      # W/m2.K
        #h_conv            = 200.0*1e-6    # W/m2.K for water
        self.hr_rad             = 5.9      # W/m2.K See example 3.12 Incropera

        # R_arm              = 0.03          # m

        #self.Tblood             = 36.6          # C blood temeprature
        self.Tair               = 33.7          # C
        self.Tinit_tissue       = 36.8 
        self.Tinit_skin         = 33.7
        self.Tinit_blood        = 37.0
        self.Tinlet_bl          = 37.0

        self.skin_L             = 0.0 # no clothing resistance if 0.

        self.w_perfusion        = 5e-4          # 1/s terminal blood flow per volume of tissue.

        self.qm_0                    = 700.0           # Basal metabolic heat
        self.Nusselt                 = 4.0#*1.0*0.42/(0.42*1.0+4*0.5) # Tt fixed rather than Tw. Nu_new=Nu* kt/(kt+Nu*kb)  

        self.tissue_input_files_set(self.tissueMeshNumber)



    def coupled_set(self, value):
        self.testFlow           = (value == 1)
        self.bioheat            = (value == 2)
        self.CoupledBioheatFlow = (value == 3)       

    def tissue_input_files_set(self, value):
        if value==1:
            self.numberOfNodesTissue    = 58557
            self.numberOfElementsTissue = 253489
            self.tissueElementsFile = 'input/bioheat/elements2.csv'
            self.tissueNodesFile    = "input/bioheat/nodes.csv"
            self.boundaryNodesFile  =  'input/bioheat/boundary_nodes.csv'
        elif value==2:
            self.numberOfNodesTissue    = 48083
            self.numberOfElementsTissue = 248166
            self.tissueElementsFile = 'input/bioheat/mesh2/elements.csv'
            self.tissueNodesFile    = "input/bioheat/mesh2/nodes.csv"
            self.boundaryNodesFile  = 'input/bioheat/mesh2/boundary_nodes.csv'
        elif value==3:
            self.numberOfNodesTissue    = 225
            self.numberOfElementsTissue = 128
            self.tissueElementsFile = 'input/bioheat/mesh3/elements.csv' 
            self.tissueNodesFile    = "input/bioheat/mesh3/nodes.csv"
            self.boundaryNodesFile  = 'input/bioheat/mesh3/boundary_nodes.csv'
        elif value==4:
            self.numberOfNodesTissue    = 10143
            self.numberOfElementsTissue = 33247
            self.tissueElementsFile = 'input/bioheat/mesh4/elements.csv' 
            self.tissueNodesFile    = "input/bioheat/mesh4/nodes.csv"
            self.boundaryNodesFile  = 'input/bioheat/mesh4/boundary_nodes.csv' 			

    def number_of_element_nodes_set(self, value): 
        if value == 3:
            self.numberOfElementNodes = 8
        else:
            self.numberOfElementNodes = 4


