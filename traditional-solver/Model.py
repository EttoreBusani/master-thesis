import numpy as np
import pyomo.environ as pyo

class PlatformModel:

    def __init__(self):
        self.eps=1e-12

#_______________________________________________________________________________
#--------------------------- PYOMO MODEL DEFINITION ----------------------------
    def create_platform_model(self, input_file):

        model = pyo.AbstractModel()
        #-----------------------------------------------------------------------
        # SETS
        #-----------------------------------------------------------------------
        model.Users = pyo.Set()
        model.Applications = pyo.Set()
        model.Users_per_Application = pyo.Set(model.Applications,within=model.Users)

        #-----------------------------------------------------------------------
        # PARAMETERS
        #-----------------------------------------------------------------------
        ########## System and applications parameters
        # c_e & c_c & N_e
        model.cloud_VM_cost=pyo.Param(default=0.0,mutable=True) # c_e
        model.edge_VM_cost=pyo.Param(default=0.0,mutable=True) # c_c
        model.max_edge_VM_number=pyo.Param(default=0.0,mutable=True)
        # D_e & D_c
        model.edge_demand_vector=pyo.Param(model.Applications,default=0.0,mutable=True) # D_e
        model.cloud_demand_vector=pyo.Param(model.Applications,default=0.0,mutable=True) # D_c
        # delta^a & lamda^a & R^a & gamma^a & r_1^a
        model.data_size=pyo.Param(model.Applications,default=0.0,mutable=True)
        model.Lambda=pyo.Param(model.Applications,default=0.0,mutable=True)
        model.R_constraints=pyo.Param(model.Applications,default=0.0,mutable=True)
        model.gamma = pyo.Param(model.Applications,default=0.0,mutable=True)
        model.r_1 = pyo.Param(model.Applications,default=0.0,mutable=True)
        # r_min & r_max & big M & T
        model.M = pyo.Param(default=10000,mutable=True)
        model.r_max = pyo.Param(default=0.0,mutable=True)
        model.r_min = pyo.Param(default=0.0,mutable=True)
        model.T = pyo.Param(default=1,mutable=True)

        ############ Users' parameters
        model.s_1=pyo.Param(model.Users, default=0.0)
        model.s_2=pyo.Param(model.Users, default=0.0)
        model.user_demand1_vector=pyo.Param(model.Users,default=0.0,mutable=True)
        model.user_demand2_vector=pyo.Param(model.Users,default=0.0,mutable=True)
        model.time=pyo.Param(model.Users,default=0.0,mutable=True)
        model.network_Bandwidth=pyo.Param(model.Users, default=0.0,mutable=True)
        model.U=pyo.Param(model.Users,default=0.0,mutable=True)
        model.beta=pyo.Param(model.Users,default=0.0,mutable=True)
        model.data_cost=pyo.Param(model.Users,default=0.0,mutable=True)
        model.alpha=pyo.Param(model.Users,default=0.0,mutable=True)
        model.p_1=pyo.Param(model.Users,default=0.0,mutable=True)
        model.p_2=pyo.Param(model.Users,default=0.0,mutable=True)

        # model.users_energy_consumption=pyo.Param(model.Users,2,default=0.0,mutable=True)
        # data transder size per application?
        # Memory? Gamma and r_1

        #-----------------------------------------------------------------------
        # VARIABLES
        #-----------------------------------------------------------------------
        #### Platform variables
        model.edge_VM_number=pyo.Var(model.Applications,within=pyo.NonNegativeIntegers)
        model.cloud_VM_number=pyo.Var(model.Applications,within=pyo.NonNegativeIntegers)
        model.y_e=pyo.Var(model.Users,within=pyo.Boolean)
        model.y_c=pyo.Var(model.Users,within=pyo.Boolean)
        model.r = pyo.Var(within=pyo.NonNegativeReals)
        #### Users' variables
        model.x_1=pyo.Var(model.Users,initialize=0,within=pyo.Binary)
        model.x_2=pyo.Var(model.Users,initialize=0,within=pyo.Binary)
        model.t_1=pyo.Var(model.Users,initialize=0,within=pyo.Binary)
        model.t_2=pyo.Var(model.Users,initialize=0,within=pyo.Binary)
        model.z=pyo.Var(model.Users,initialize=0,within=pyo.Binary)   # z_i^{(12)}

        #-----------------------------------------------------------------------
        # CONSTRAINTS
        #-----------------------------------------------------------------------
        # Constraint on max Edge nodes available
        model.edge_VM_number_constraints = pyo.Constraint(rule=self.max_edge_VM_number_constraint_MEC)
        # model.cloud_VM_number_constraint = pyo.Constraint(rule=self.max_cloud_VM_number_constraint_MEC)

        # Edge and cloud load constraints
        model.edge_saturation_constraints = pyo.Constraint(model.Applications, rule=self.edge_saturation_constraint_MEC)
        model.cloud_saturation_constraints = pyo.Constraint(model.Applications, rule=self.cloud_saturation_constraint_MEC)

        # Response time constraints
        model.response_time_constraints = pyo.Constraint(model.Applications,model.Users,rule=self.response_time_constraint_MEC1)

        # User's assignment constraints
        model.y_constraints = pyo.Constraint(model.Users,rule=self.y_constraint_MEC)

        # Linear constraints for Users'problem
        model.linear_constraint_t1_lb = pyo.Constraint(model.Applications,model.Users,rule=self.linear_constraint_t1_lb)
        model.linear_constraint_t1_ub = pyo.Constraint(model.Applications,model.Users,rule=self.linear_constraint_t1_ub)
        model.linear_constraint_t2_lb = pyo.Constraint(model.Applications,model.Users,rule=self.linear_constraint_t2_lb)
        model.linear_constraint_t2_ub = pyo.Constraint(model.Applications,model.Users,rule=self.linear_constraint_t2_ub)
        model.linear_constraint_z_lb = pyo.Constraint(model.Applications,model.Users,rule=self.linear_constraint_z_lb)
        model.linear_constraint_z_ub = pyo.Constraint(model.Applications,model.Users,rule=self.linear_constraint_z_ub)

        # model.linear_constraint_t1 = pyo.Constraint(model.Applications,model.Users,rule=self.alternative_t1)
        # model.linear_constraint_t2 = pyo.Constraint(model.Applications,model.Users,rule=self.alternative_t2)
        # model.linear_constraint_z = pyo.Constraint(model.Applications,model.Users,rule=self.alternative_z)


        model.linear_constraint_x1_lb = pyo.Constraint(model.Users,rule=self.linear_constraint_x1_lb)
        model.linear_constraint_x1_ub = pyo.Constraint(model.Users,rule=self.linear_constraint_x1_ub)
        model.linear_constraint_x2_lb = pyo.Constraint(model.Users,rule=self.linear_constraint_x2_lb)
        model.linear_constraint_x2_ub = pyo.Constraint(model.Users,rule=self.linear_constraint_x2_ub)
        model.one_depl_constraint = pyo.Constraint(model.Users,rule=self.one_depl_constraint)

        # Constraint on r
        model.r_constraint_min = pyo.Constraint(rule=self.r_constraint_min)
        model.r_constraint_max = pyo.Constraint(rule=self.r_constraint_max)

        #-----------------------------------------------------------------------
        # OBJECTIVE FUNCTION
        #-----------------------------------------------------------------------
        model.Maximum_edge_profit=pyo.Objective(rule=self.objective_function_MEC, sense=pyo.maximize)

        ### CALL TO LOAD PLATFORM DATA -> FROM ABSTRACT TO CONCRETE
        data=self.load_platform_data(input_file, model)

        return data, model
#-------------------------------------------------------------------------------------------------------------------
# DEFINITIONS

    def compute_cost(model,a,i,depl):
        if(depl==1):
            return model.alpha[i]*model.r_1[a]*model.time[i]/3600 + (1-model.alpha[i])*model.beta[i]/3600000*model.p_1[i]*model.Lambda[a]*model.time[i]*model.time[i]
        else:
            return model.alpha[i]*(model.r_1[a]+model.r)*model.time[i]/3600 + (1-model.alpha[i])*model.beta[i]/3600000*model.p_2[i]*model.Lambda[a]*model.time[i]*model.time[i] + model.time[i]*model.data_size[a]*model.Lambda[a]*model.data_cost[i]


    ### Objective Function Definition
    def objective_function_MEC(self, model):
        # return sum(model.time[i]/3600 * (model.x_1[i]*model.r_1[a] + model.x_2[i]*(model.r_1[a]+model.gamma[a]*model.r)) for a in model.Applications for i in model.Users_per_Application[a])-sum((model.edge_VM_cost*model.edge_VM_number[a] + model.cloud_VM_cost*model.cloud_VM_number[a])*model.T/3600 for a in model.Applications)
        revenue = 0
        cost = 0
        for a in model.Applications:
            cost += (model.edge_VM_cost*model.edge_VM_number[a] + model.cloud_VM_cost*model.cloud_VM_number[a])*model.T/3600
            for i in model.Users_per_Application[a]:
                revenue += model.time[i]/3600 * (model.x_1[i]*model.r_1[a] + model.x_2[i]*(model.r_1[a]+model.gamma[a]*model.r))
        return revenue-cost

    ### Contraints Definition
    def  max_edge_VM_number_constraint_MEC(self, model):
        return sum(model.edge_VM_number[a] for a in model.Applications)<=model.max_edge_VM_number

    def max_cloud_VM_number_constraint_MEC(self,model):
        return sum(model.cloud_VM_number[a] for a in model.Applications)<=100

    def  y_constraint_MEC(self, model, i):
        return  model.y_e[i]+model.y_c[i]==model.x_2[i]

    def  edge_saturation_constraint_MEC(self, model, a):
        return model.edge_VM_number[a]*(1-self.eps) - sum(model.edge_demand_vector[a] * model.Lambda[a] *  model.y_e[i] for i in model.Users_per_Application[a]) >= 0

    def cloud_saturation_constraint_MEC(self, model, a):
        return model.cloud_VM_number[a]*(1-self.eps) - sum(model.cloud_demand_vector[a] * model.Lambda[a] *  model.y_c[i] for i in model.Users_per_Application[a]) >= 0

    def response_time_constraint_MEC(self, model, a):
        respected = True
        L_e = model.edge_demand_vector[a]*model.Lambda[a]*sum(model.y_e[i] for i in model.Users_per_Application[a])
        L_c = model.cloud_demand_vector[a]*model.Lambda[a]*sum(model.y_c[i] for i in model.Users_per_Application[a])
        for i in model.Users_per_Application[a]:
            term1 = model.user_demand2_vector[i]*model.x_2[i]
            term2 = model.data_size[a]*model.x_2[i]/model.network_Bandwidth[i]
            term3 = model.edge_demand_vector[a]*model.y_e[i]/(1-L_e/(model.edge_VM_number[a]+self.eps)) + model.cloud_demand_vector[a]*model.y_c[i]/(1-L_c/(model.cloud_VM_number[a]+self.eps))
            local_respected = term1+term2+term3 <= model.R_constraints[a]
            respected = respected and local_respected
        return respected

    def response_time_constraint_MEC1(self,model,a,i):
        if(i in model.Users_per_Application[a]):
            L_e = model.edge_demand_vector[a]*model.Lambda[a]*sum(model.y_e[i] for i in model.Users_per_Application[a])
            L_c = model.cloud_demand_vector[a]*model.Lambda[a]*sum(model.y_c[i] for i in model.Users_per_Application[a])
            term1 = model.user_demand2_vector[i]*model.x_2[i]
            term2 = model.data_size[a]*model.x_2[i]/model.network_Bandwidth[i]
            term3 = model.edge_demand_vector[a]*model.y_e[i]/(1-L_e/(model.edge_VM_number[a]+self.eps)) + model.cloud_demand_vector[a]*model.y_c[i]/(1-L_c/(model.cloud_VM_number[a]+self.eps))
            return term1+term2+term3 <= model.R_constraints[a]
        else: return pyo.Constraint.Skip
        # return sum(model.user_demand_matrix[i,k] * model.x[i,k] for k in model.Deployments)+\
        #     sum(model.data_size[k] * model.x[i,k]/model.network_Bandwidth[i] for k in model.Deployments if k != "s1")+\
        #     sum((model.edge_VM_number[k] * model.edge_demand_matrix[k] * model.x[i,k] * model.y_e[i])/(model.edge_VM_number[k]-sum(model.edge_demand_matrix[k] * model.x[i1,k] * model.Lambda *  model.y_e[i1] for i1 in model.Users) + self.eps)for k in model.Deployments if k != "s1")+\
        #         sum((model.cloud_VM_number[k] * model.cloud_demand_matrix[k] * model.x[i,k] * model.y_c[i])/(model.cloud_VM_number[k]-sum(model.cloud_demand_matrix[k] * model.x[i1,k] * model.Lambda *  model.y_c[i1] for i1 in model.Users) + self.eps)for k in model.Deployments if k != "s1") <= model.R_constraints

    def r_constraint_min(self, model):
        return model.r>=model.r_min

    def r_constraint_max(self, model):
        return model.r<=model.r_max

    ## Linear Constraints to replace Users' optimization problem ##
    def alternative_t1(self,model,a,i):
        if(i in model.Users_per_Application[a]):
            cost = model.alpha[i]*model.r_1[a]*model.time[i]/3600 + (1-model.alpha[i])*model.beta[i]/3600000*model.p_1[i]*model.Lambda[a]*model.time[i]*model.time[i]
            return model.t_1[i] == (model.U[i]/3600*model.time[i]-cost >= 0)
        else: return pyo.Constraint.Skip

    def alternative_t2(self,model,a,i):
        if(i in model.Users_per_Application[a]):
            cost = model.alpha[i]*(model.r_1[a]+model.r)*model.time[i]/3600 + (1-model.alpha[i])*model.beta[i]/3600000*model.p_2[i]*model.Lambda[a]*model.time[i]*model.time[i] + model.time[i]*model.data_size[a]*model.Lambda[a]*model.data_cost[i]
            return model.t_2[i] == (model.U[i]/3600*model.time[i]-cost >= 0)
        else: return pyo.Constraint.Skip

    def alternative_z(self,model,a,i):
        if(i in model.Users_per_Application[a]):
            cost1 = model.alpha[i]*model.r_1[a]*model.time[i]/3600 + (1-model.alpha[i])*model.beta[i]/3600000*model.p_1[i]*model.Lambda[a]*model.time[i]*model.time[i]
            cost2 = model.alpha[i]*(model.r_1[a]+model.r)*model.time[i]/3600 + (1-model.alpha[i])*model.beta[i]/3600000*model.p_2[i]*model.Lambda[a]*model.time[i]*model.time[i] + model.time[i]*model.data_size[a]*model.Lambda[a]*model.data_cost[i]
            return model.z[i] == (cost2 - cost1 >= 0)
        else: return pyo.Constraint.Skip

    def  linear_constraint_t1_lb(self, model, a, i):
        if(not i in model.Users_per_Application[a]): return pyo.Constraint.Skip
        else:
            cost = model.alpha[i]*model.r_1[a]*model.time[i]/3600 + (1-model.alpha[i])*model.beta[i]/3600000*model.p_1[i]*model.Lambda[a]*model.time[i]*model.time[i]
            return (-1)*(1-model.t_1[i])*model.M <= model.U[i]/3600*model.time[i]-cost

    def  linear_constraint_t1_ub(self, model, a, i):
        if(not i in model.Users_per_Application[a]): return pyo.Constraint.Skip
        else:
            cost = model.alpha[i]*model.r_1[a]*model.time[i]/3600 + (1-model.alpha[i])*model.beta[i]/3600000*model.p_1[i]*model.Lambda[a]*model.time[i]*model.time[i]
            return model.U[i]/3600*model.time[i]-cost <= model.t_1[i]*model.M

    def  linear_constraint_t2_lb(self, model, a, i):
        if(not i in model.Users_per_Application[a]): return pyo.Constraint.Skip
        else:
            cost = model.alpha[i]*(model.r_1[a]+model.r)*model.time[i]/3600 + (1-model.alpha[i])*model.beta[i]/3600000*model.p_2[i]*model.Lambda[a]*model.time[i]*model.time[i] + model.time[i]*model.data_size[a]*model.Lambda[a]*model.data_cost[i]
            return (-1)*(1-model.t_2[i])*model.M <= model.U[i]/3600*model.time[i]-cost

    def  linear_constraint_t2_ub(self, model, a, i):
        if(not i in model.Users_per_Application[a]): return pyo.Constraint.Skip
        else:
            cost = model.alpha[i]*(model.r_1[a]+model.r)*model.time[i]/3600 + (1-model.alpha[i])*model.beta[i]/3600000*model.p_2[i]*model.Lambda[a]*model.time[i]*model.time[i] + model.time[i]*model.data_size[a]*model.Lambda[a]*model.data_cost[i]
            return model.U[i]/3600*model.time[i]-cost <= model.t_2[i]*model.M

    # z_i is equal to one if it better for user i to run depl 1 (i.e. cost1 < cost2)
    def linear_constraint_z_lb(self, model, a, i):
        if(not i in model.Users_per_Application[a]): return pyo.Constraint.Skip
        else:
            cost1 = model.alpha[i]*model.r_1[a]*model.time[i]/3600 + (1-model.alpha[i])*model.beta[i]/3600000*model.p_1[i]*model.Lambda[a]*model.time[i]*model.time[i]
            cost2 = model.alpha[i]*(model.r_1[a]+model.r)*model.time[i]/3600 + (1-model.alpha[i])*model.beta[i]/3600000*model.p_2[i]*model.Lambda[a]*model.time[i]*model.time[i] + model.time[i]*model.data_size[a]*model.Lambda[a]*model.data_cost[i]
            return (-1)*(1-model.z[i])*model.M <= cost2-cost1

    def linear_constraint_z_ub(self, model, a, i):
        if(not i in model.Users_per_Application[a]): return pyo.Constraint.Skip
        else:
            cost1 = model.alpha[i]*model.r_1[a]*model.time[i]/3600 + (1-model.alpha[i])*model.beta[i]/3600000*model.p_1[i]*model.Lambda[a]*model.time[i]*model.time[i]
            cost2 = model.alpha[i]*(model.r_1[a]+model.r)*model.time[i]/3600 + (1-model.alpha[i])*model.beta[i]/3600000*model.p_2[i]*model.Lambda[a]*model.time[i]*model.time[i] + model.time[i]*model.data_size[a]*model.Lambda[a]*model.data_cost[i]
            return cost2-cost1 <= model.z[i]*model.M

    def linear_constraint_x1_ub(self, model, i):
        return model.x_1[i] <= model.s_1[i]*model.t_1[i]

    def linear_constraint_x1_lb(self, model, i):
        return model.x_1[i] >= model.s_1[i]*model.t_1[i] + model.s_2[i]*(model.z[i]-1)

    def linear_constraint_x2_ub(self, model, i):
         return model.x_2[i] <= model.s_2[i]*model.t_2[i]

    def linear_constraint_x2_lb(self, model, i):
         return model.x_2[i] >= model.s_2[i]*model.t_2[i] - model.s_1[i]*model.z[i]

    def one_depl_constraint(self, model, i):
        return model.x_1[i]+model.x_2[i] <= 1


############################### End Pyomo Model ########################################


    ## Function to load the information of model from .dat file
    def load_platform_data(self,file_name,PlatformModel):
        ### Data Object
        PlatformData = pyo.DataPortal(model=PlatformModel)

        ### SETS
        PlatformData.load(filename=file_name, Set=PlatformModel.Users, format="set")
        PlatformData.load(filename=file_name, Set=PlatformModel.Applications, format="set")
        PlatformData.load(filename=file_name, Set=PlatformModel.Users_per_Application, format="set_array")

        ### SYSTEM AND APPLICATION PARAMETERS
        PlatformData.load(filename=file_name, Param=PlatformModel.cloud_VM_cost, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.edge_VM_cost, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.max_edge_VM_number, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.edge_demand_vector, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.cloud_demand_vector, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.data_size, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.Lambda, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.R_constraints, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.gamma, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.r_1, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.M, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.r_max, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.r_min, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.T, format="param")

        ### USERS' parameters
        PlatformData.load(filename=file_name, Param=PlatformModel.s_1, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.s_2, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.user_demand1_vector, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.user_demand2_vector, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.time, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.network_Bandwidth, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.U, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.beta, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.data_cost, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.alpha, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.p_1, format="param")
        PlatformData.load(filename=file_name, Param=PlatformModel.p_2, format="param")

        return PlatformData


# class UserModel:
#
#     def __init__(self, S, alpha_treatment):
#
#         self.generate_random_user_parameters(S, alpha_treatment)
#         #self.min_cost = np.inf
#         #self.max_cost = 0
#
#
#     ## Function to generate the users' data randomly
#     def generate_random_user_parameters(self, S, alpha_treatment):
#
#         kwh_to_ws=1/3600000
#         beta=np.array(np.random.uniform(0.491*kwh_to_ws, 0.521*kwh_to_ws,S.N),dtype=float)     # user's beta
#         zeta=np.array(np.random.uniform(5.3e-5, 6e-5, S.N),dtype=float)
#         U=np.array([(1+(i >= S.N/2))*np.random.uniform(2,2.5) for i in range(S.N)],dtype=float)
#
#         # import pdb; pdb.set_trace()
#
#         E=np.array([0 for i in range(S.N)],dtype=float)     # energy
#         M=np.array([0 for i in range(S.N)],dtype=int)   # memory
#         e=np.array([[0 for it2 in range(S.D)] for it1 in range(S.N)],dtype=float)   # matrix of p_i^{(k)}
#         demands=np.array([[0 for it2 in range(S.D)] for it1 in range(S.N)],dtype=float)   # matrix of D_i^{(k)}
#         T_i=np.array([np.random.uniform(540,660)*(1+(i >= S.N/2)) for i in range(S.N)],dtype=float)
#         B=np.array([0 for i in range(S.N)],dtype=float)
#
#         demands_list = [[3.211, 0.52171], [3.211, 0.7424, 0.19131], [3.211, 0.7424, 0.52171, 0.19131], [3.211, 1.1054, 0.7424, 0.52171, 0.19131]]
#
#         for i in range(S.N):
#             for k in range(0,S.D):
#                 demands[i][k]=round(np.random.uniform(demands_list[S.D-2][k]*0.95, demands_list[S.D-2][k]*1.05),5)
#                 if k == 0:
#                     e[i][k]=demands[i][k]*5
#                 else:
#                     e[i][k]=demands[i][k]*10
#
#             B[i]=round(np.random.uniform(9.9,25.1),2)/8 # MB/s
#             lower_bound_E = S.Lambda*T_i[i]*T_i[i]*e[i][-1]*0.9
#             upper_bound_E = S.Lambda*T_i[i]*T_i[i]*e[i][0]*1.2
#             E[i]=round(np.random.uniform(lower_bound_E, upper_bound_E),3)
#             M[i]=np.random.randint(S.D*5,S.D*6)
#
#         if alpha_treatment == 'random':
#             alpha=np.array(np.random.uniform(0.5,0.75,S.N),dtype=float)
#         elif alpha_treatment == 'fixed':
#             alpha=np.array([0.65 for jj in range(S.N)], dtype=float)
#         elif alpha_treatment == 'energy_dependent':
#             alpha=np.array([0 for i in range(S.N)],dtype=float)
#             for i in range(S.N):
#                 alpha[i] = 0.3*(E[i] - lower_bound_E)/(upper_bound_E-lower_bound_E) + 0.5
#         else:
#             print("Wrong alpha treatment")
#
#         self.beta=beta
#         self.zeta=zeta
#         self.alpha=alpha
#         self.U = U
#         self.max_energy=E
#         self.max_memory=M
#         self.dep_power_consumption=e
#         self.demands=demands
#         self.B=B
#         self.T=T_i
#
#
#     ## Function to check if the energy and memory of users' device is enough to run the task for users
#     def users_memory_energy_checking(self, S):
#
#         s=np.array([[0 for it2 in range(S.D)] for it1 in range(S.N)],dtype=float)
#
#         for i in range(S.N):
#             for k in range(S.D):
#                 if S.Lambda*self.T[i]*self.T[i]*self.dep_power_consumption[i][k]<=self.max_energy[i] and S.memory[k]<=self.max_memory[i] :
#                     s[i][k]=1
#
#         return s
#
#
#     ## Function to solve the users' problem
#     def solve_users_problem(self, S, r):
#
#         best_solution=np.array([-1 for i in range(S.N)],dtype=int)
#         x=np.array([[0 for it2 in range(S.D)] for it1 in range(S.N)],dtype=float)
#
#         for i in range(S.N):
#             best_cost = self.U[i]
#
#             for k in range(S.D):
#
#                 cost = self.T[i]*(self.alpha[i]*r[k] + (1-self.alpha[i])*self.beta[i]*self.dep_power_consumption[i][k]*S.Lambda*self.T[i] + self.zeta[i]*S.data_size[k]*S.Lambda)
#
#                 #if cost < self.min_cost:
#                 #    self.min_cost = cost
#                 #if cost > self.max_cost:
#                 #    self.max_cost = cost
#                 #print(self.T[i]*(self.beta[i]*self.dep_power_consumption[i][k]*S.Lambda*self.T[i] + S.zeta*S.data_size[k]*S.Lambda))
#                 #print(cost)
#
#                 #print("r1 = " + str(r[0]))
#                 #print("r_param = " + str(r[1]-r[0]))
#                 #print("r = " + str(r))
#                 #print(self.alpha[i]*r[k])
#                 #print((1-self.alpha[i])*self.beta[i]*self.dep_power_consumption[i][k]*S.Lambda*self.T[i])
#                 #print(self.zeta[i]*S.data_size[k]*S.Lambda)
#                 #print("\n")
#
#                 if cost <= best_cost and (S.Lambda*self.T[i]*self.T[i]*self.dep_power_consumption[i][k]<=self.max_energy[i] and S.memory[k]<=self.max_memory[i]):
#                     best_cost=cost
#                     best_solution[i]=k
#
#         for i in range(S.N):
#            if best_solution[i]>-1:
#                x[i][best_solution[i]]=1
#
#         return x
