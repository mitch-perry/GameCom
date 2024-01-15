import numpy as np
import scipy.linalg as linalg
import scipy.sparse as sparse
import scipy.sparse.linalg as sp_linalg
import cvxpy as cp
import gurobipy
import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter('ignore')

# Set solver to use with cvxpy. See "Choosing a solver" section here for using alternatives 
# to Gurobi: https://www.cvxpy.org/tutorial/advanced/index.html
cp_solver = 'GUROBI'

# Functions provided here are generalizations of the 
# functions defined in the two notebooks "ecoli_model_compute_ss_gne.ipynb"
# and "ecoli_model_stability.ipynb" found in the folder
# EColiModel/GameCom/.

def initial_guess(target_bm, species_params, S, Jl, lumen_reactions_idx, 
                  death_rate, reaction_lb, reaction_ub):
    solutions = [cp.Variable((species_params[i]['J'] + Jl,1)) for i in range(len(species_params))]
    constraints = []
    ex_constraint = 0
    for i in range(len(species_params)):
        constraints.append(S[:,np.concatenate([species_params[i]['reactions_idx'], lumen_reactions_idx]).flatten()] @ solutions[i] == 0)
        constraints.append(death_rate <= species_params[i]['e'].T @ solutions[i])
        constraints.append(reaction_lb[species_params[i]['reactions_idx'].flatten()] <= solutions[i][0:species_params[i]['J']])
        constraints.append(reaction_ub[species_params[i]['reactions_idx'].flatten()] >= solutions[i][0:species_params[i]['J']])
        ex_constraint = ex_constraint + target_bm[i] * solutions[i][species_params[i]['J']:]
    constraints.append(ex_constraint >= reaction_lb[lumen_reactions_idx.flatten()])
    constraints.append(ex_constraint <= reaction_ub[lumen_reactions_idx.flatten()])
    prob = cp.Problem(cp.Maximize(0), constraints)
    prob.solve(solver=cp_solver)
    return [prob.status] + [solutions[i].value for i in range(len(species_params))]

def compute_steady_state(target_bm, x_values, max_iters, tol_bm, tol_flux,
                        delta_vals, species_params, S, Jl, lumen_reactions_idx,
                        reaction_lb, reaction_ub, death_rate):
    iteration = 0
    bm_vals = target_bm
    x_vals = [x_values[i] * bm_vals[i] for i in range(len(bm_vals))] # These are biomass-weighted fluxes.
    max_current_change_bm = 1
    max_current_change = 1

    #while iteration < max_iters and (max_current_change_bm > tol_bm or max_current_change > tol_flux):
    while iteration < max_iters and max_current_change_bm > tol_bm:
        print('Iteration ', iteration)
        delta = delta_vals[iteration]
        solution_flux_vals = []
        solution_bm_vals = []

        for i in range(len(target_bm)):
            solution = cp.Variable((species_params[i]['J'] + Jl + 1, 1))
            constraints = []
            constraints.append(S[:, np.concatenate([species_params[i]['reactions_idx'], lumen_reactions_idx]).flatten()] @ solution[0:-1] == 0)
            constraints.append(solution[0:species_params[i]['J']] >= solution[-1] * reaction_lb[species_params[i]['reactions_idx'].flatten()])
            constraints.append(solution[0:species_params[i]['J']] <= solution[-1] * reaction_ub[species_params[i]['reactions_idx'].flatten()])
            constraints.append(solution[-1] * death_rate <= species_params[i]['e'].T.dot(solution[0:-1]))
            ex_constraint = 0
            for j in range(len(target_bm)):
                if j != i:
                    ex_constraint = ex_constraint + x_vals[j][species_params[j]['J']:]
                else:
                    ex_constraint = ex_constraint + solution[species_params[i]['J']:-1]
            constraints.append(ex_constraint <= reaction_ub[lumen_reactions_idx.flatten()])
            constraints.append(ex_constraint >= reaction_lb[lumen_reactions_idx.flatten()])
            objective = solution[-1] - 0.5 * delta * cp.quad_form(solution[0:species_params[i]['J']+Jl] - x_vals[i], sparse.identity((species_params[i]['J']+Jl)).tocsr()) - 0.5 * delta * cp.power(solution[-1] - bm_vals[i], 2)
            prob = cp.Problem(cp.Maximize(objective), constraints)
            num_tries = 0
            while prob.status is None and num_tries < 2:
                try:
                    prob.solve(solver=cp_solver)
                    num_tries = num_tries + 1
                except:
                    num_tries = num_tries + 1
            if solution.value is None:
                #print('solution.value was none')
                solution.value = np.zeros((species_params[i]['J']+Jl+1,1))
            if prob.status != 'optimal':
                #print('problem.status was opt-infeasible')
                solution.value = np.zeros((species_params[i]['J']+Jl+1,1))
            if abs(solution.value[-1] - bm_vals[i]) / np.sum(bm_vals) > tol_bm:
                solution_bm_vals.append(solution.value[-1])
                solution_flux_vals.append(solution.value[0:-1])
            else:
                solution_bm_vals.append(bm_vals[i])
                solution_flux_vals.append(x_vals[i])
        
        # Update fluxes and biomasses.
        x_old_vals = np.copy(x_vals)
        bm_old_vals = np.copy(bm_vals)
        for i in range(len(bm_vals)):
            x_vals[i] = x_vals[i] + (1 / (10 + np.sqrt(iteration))) * (solution_flux_vals[i] - x_vals[i])
            bm_vals[i] = bm_vals[i] + (1 / (2 + np.sqrt(iteration))) * (solution_bm_vals[i] - bm_vals[i])
        
        # Re-scale biomass values so community has total biomass of 1.
        bm_vals = np.array(bm_vals) / np.sum(bm_vals)

        max_current_change = 0
        max_current_change_bm = 0
        for i in range(len(target_bm)):
            if np.linalg.norm(x_old_vals[i]) > 0:
                change_i = np.linalg.norm(x_vals[i] - x_old_vals[i]) / np.linalg.norm(x_old_vals[i])
                if change_i > max_current_change:
                    max_current_change = change_i
                
            if bm_old_vals[i] > 0:
                change_i_bm = abs(bm_vals[i] - bm_old_vals[i]) / np.sum(bm_old_vals)
                if change_i_bm > max_current_change_bm:
                    max_current_change_bm = change_i_bm
        
        iteration = iteration + 1
        print('max_current_change_bm: ', max_current_change_bm)
        #print('max_current_change: ', max_current_change)
        print('bm_vals: ', bm_vals)
    
    # Convert fluxes back to being non-biomass-weighted.
    x_vals_new = []
    for i in range(len(bm_vals)):
        if bm_vals[i] == 0:
            x_vals_new.append(np.zeros(x_vals[i].shape))
        else:
            x_vals_new.append(x_vals[i] / bm_vals[i])
    
    return [iteration, x_vals_new, bm_vals]


def stability(bm_vals, x_vals, pert_size, species_params, S, Jl, 
             reaction_lb, reaction_ub, lumen_reactions_idx, 
             lumen_metabolites_idx, feas_tol):
    print('Starting stability function')
    # Re-run best response problem with Gurobi/simplex method to get optimal basis.
    U_vals = []
    L_vals = []
    F_vals = []
    lambdas = []
    lambdas_L = []
    lambdas_U = []
    rxns_vals = []
    rxns_ex_vals = []

    for i in range(len(bm_vals)):
        #print('Re-solving problem for species ', i)
        # Species i's problem.
        if bm_vals[i] == 0:
            U_vals.append(np.array([]))
            L_vals.append(np.array([]))
            F_vals.append(np.array([]))
            lambdas.append(np.zeros(Jl,))
            lambdas_L.append(np.zeros((species_params[i]['J'],)))
            lambdas_U.append(np.zeros((species_params[i]['J'],)))
            rxns_vals.append(np.zeros((species_params[i]['J'],)))
            rxns_ex_vals.append(np.zeros((Jl,)))

        else:
            env = gurobipy.Env(empty=True)
            env.setParam("OutputFlag",0)
            env.start()
            m = gurobipy.Model(species_params[i]['name'], env=env)
            
            # Create variables.
            lb_i = reaction_lb[species_params[i]['reactions_idx'].flatten()]
            ub_i = reaction_ub[species_params[i]['reactions_idx'].flatten()]

            lb_i_pert = np.random.random((species_params[i]['J'], 1))
            lb_il_pert = np.random.random((Jl, 1))
            ub_i_pert = np.random.random((species_params[i]['J'], 1))

            lb_i = lb_i - pert_size * lb_i_pert
            ub_i = ub_i + pert_size * ub_i_pert

            rxns_i = m.addMVar(shape=species_params[i]['J'], name='rxns'+str(i), lb=lb_i.flatten(), ub=ub_i.flatten())

            lb_i_ex = reaction_lb[lumen_reactions_idx.flatten()]
            ub_i_ex = reaction_ub[lumen_reactions_idx.flatten()]

            for j in range(len(bm_vals)):
                if j != i:
                    lb_i_ex = lb_i_ex - bm_vals[j] * x_vals[j][species_params[j]['J']:]
                    ub_i_ex = ub_i_ex - bm_vals[j] * x_vals[j][species_params[j]['J']:]
            
            if bm_vals[i] != 0:
                lb_i_ex = lb_i_ex / bm_vals[i]
                ub_i_ex = ub_i_ex / bm_vals[i]
            
            lb_i_ex = lb_i_ex - pert_size * lb_il_pert

            rxns_i_ex = m.addMVar(shape=Jl, name='rxns_ex_'+str(i), lb=lb_i_ex.flatten(), ub=ub_i_ex.flatten());

            # Set objective.
            m.setObjective(species_params[i]['e'].toarray()[0:species_params[i]['J']].T @ rxns_i, gurobipy.GRB.MAXIMIZE);

            # Make constraints.
            S_i = S[np.concatenate([species_params[i]['metabolites_idx'], lumen_metabolites_idx]).flatten(), :]
            S_i_ex = S_i[:, lumen_reactions_idx.flatten()]
            S_i = S_i[:, species_params[i]['reactions_idx'].flatten()]
            m.addConstr(S_i @ rxns_i + S_i_ex @ rxns_i_ex == np.zeros((S_i.shape[0],)), name='internal_fba'+str(i));

            # Solve.
            m.params.Method = 0;
            m.params.FeasibilityTol = feas_tol;
            m.update();
            m.optimize();

            #print('m.status: ', m.status)

            # Check degeneracy.
            primal_basis_i = rxns_i.VBasis == 0
            primal_basis_i_ex = rxns_i_ex.VBasis == 0

            active_basic_i_lb = np.where(rxns_i.X[primal_basis_i] - lb_i[primal_basis_i] == 0)[0]
            active_basic_i_ub = np.where(ub_i[primal_basis_i] - rxns_i.X[primal_basis_i] == 0)[0]
            active_basic_i_ex_lb = np.where(rxns_i_ex.X[primal_basis_i_ex] - lb_i_ex[primal_basis_i_ex] == 0)[0]
            active_basic_i_ex_ub = np.where(ub_i_ex[primal_basis_i_ex] - rxns_i_ex.X[primal_basis_i_ex] == 0)[0]

            if len(active_basic_i_lb) == 0 and len(active_basic_i_ub) == 0 and len(active_basic_i_ex_lb) == 0 and len(active_basic_i_ex_ub) == 0:
                degenerate = False
            else:
                degenerate = True

            if degenerate:
                raise Exception('Problem ' + str(i) + ' is degenerate')
            
            U_vals.append(np.where(ub_i - rxns_i.X < 1e-6)[0])
            L_vals.append(np.where(rxns_i.X - lb_i < 1e-6)[0])
            F_vals.append(np.where(np.logical_and(~np.isin(np.arange(species_params[i]['J']), U_vals[i]),
                                                ~np.isin(np.arange(species_params[i]['J']), L_vals[i])))[0])
            lambdas.append(rxns_i_ex.RC)

            lambdas_L_i = np.zeros((len(rxns_i.X),))
            lambdas_L_i[np.where(rxns_i.VBasis == -1)[0]] = rxns_i.RC[np.where(rxns_i.VBasis == -1)[0]]
            lambdas_L.append(lambdas_L_i)

            lambdas_U_i = np.zeros((len(rxns_i.X),))
            lambdas_U_i[np.where(rxns_i.VBasis == -2)[0]] = rxns_i.RC[np.where(rxns_i.VBasis == -2)[0]]
            lambdas_U.append(lambdas_U_i)

            rxns_vals.append(rxns_i.X)
            rxns_ex_vals.append(rxns_i_ex.X)
    
    E = -reaction_lb[lumen_reactions_idx.flatten()].flatten()
    for i in range(len(bm_vals)):
        E = E + bm_vals[i] * rxns_ex_vals[i]
    E = np.where(E < 1e-6)[0]

    for i in range(len(bm_vals)):
        x_vals[i][0:species_params[i]['J']] = rxns_vals[i][:,None]
        x_vals[i][species_params[i]['J']:] = rxns_ex_vals[i][:,None]
    
    b_vals = []
    for i in range(len(bm_vals)):
        b = reaction_lb[lumen_reactions_idx.flatten()]
        for j in range(len(bm_vals)):
            if j != i:
                b = b - bm_vals[j] * x_vals[j][species_params[j]['J']:]
        b = b / bm_vals[i]
        b_vals.append(b)
    
    num_vars = len(species_params)**2 + (len(species_params)**2)*Jl + len(species_params)*np.sum([F_vals[l].shape[0] for l in range(len(species_params))])
    num_constraints = len(species_params)**2 + (len(species_params)**2)*E.shape[0] + len(species_params)*np.sum([len(species_params[l]['metabolites_idx']) for l in range(len(species_params))]) + (len(species_params)**2)*Jl

    A_stable = sparse.lil_matrix((num_constraints, num_vars))
    b_stable = np.zeros((num_constraints,))

    A_stable = sparse.lil_matrix((num_constraints, num_vars))
    b_stable = np.zeros((num_constraints,))

    print('First block of constructing A')
    # Make constraints for partials of h_{k}(x) values.
    for k in range(len(species_params)):
        for j in range(len(species_params)):
            row_kj = np.zeros((num_vars,))
            row_kj[len(species_params)*k + j] = 1
            for l in range(len(species_params)):
                if l != k:
                    if bm_vals[k] > 0:
                        row_kj[(len(species_params)**2 + (len(species_params)*l + j)*Jl):(len(species_params)**2 + (len(species_params)*l + j + 1)*Jl)] = (bm_vals[l] / bm_vals[k]) * lambdas[k].T
                    else:
                        row_kj[(len(species_params)**2 + (len(species_params)*l + j)*Jl):(len(species_params)**2 + (len(species_params)*l + j + 1)*Jl)] = np.zeros(lambdas[k].T.shape)

            A_stable[len(species_params)*k + j, :] = row_kj

            if j == k:
                if bm_vals[k] > 0:
                    b_stable[len(species_params)*k + j] = -(1/bm_vals[k]) * lambdas[k].T.dot(b_vals[k].flatten())
                else:
                    b_stable[len(species_params)*k + j] = 0
            else:
                if bm_vals[k] > 0:
                    b_stable[len(species_params)*k + j] = -(1/bm_vals[k]) * lambdas[k].T.dot(x_vals[j][species_params[j]['J']:])
                else:
                    b_stable[len(species_params)*k + j] = 0

    print('Second block of constructing A')
    # Make constraints for partials of R_{k}^{E} \nu_{k}^{ex} with respect to x_{j}.
    for k in range(len(species_params)):
        for j in range(len(species_params)):
            block_kj = np.zeros((E.shape[0], num_vars))
            block_kj[:, len(species_params)**2 + (len(species_params)*k + j)*Jl + E] = np.eye(E.shape[0])

            for l in range(len(species_params)):
                if l != k:
                    if bm_vals[k] > 0:
                        block_kj[:, len(species_params)**2 + (len(species_params)*l + j)*Jl + E] = (bm_vals[l] / bm_vals[k]) * np.eye(E.shape[0])
                    else:
                        block_kj[:, len(species_params)**2 + (len(species_params)*l + j)*Jl + E] = np.zeros(E.shape[0])
            
            A_stable[(len(species_params)**2 + (len(species_params)*k + j)*E.shape[0]):(len(species_params)**2 + (len(species_params)*k + j + 1)*E.shape[0]), :] = sparse.csr_matrix(block_kj)

            if j == k:
                if bm_vals[k] > 0:
                    b_stable[(len(species_params)**2 + (len(species_params)*k + j)*E.shape[0]):(len(species_params)**2 + (len(species_params)*k + j + 1)*E.shape[0])] = -(1 / bm_vals[k]) * b_vals[k].flatten()[E]
                else:
                    b_stable[(len(species_params)**2 + (len(species_params)*k + j)*E.shape[0]):(len(species_params)**2 + (len(species_params)*k + j + 1)*E.shape[0])] = 0
            else:
                if bm_vals[k] > 0:
                    b_stable[(len(species_params)**2 + (len(species_params)*k + j)*E.shape[0]):(len(species_params)**2 + (len(species_params)*k + j + 1)*E.shape[0])] = -(1 / bm_vals[k]) * x_vals[j][species_params[j]['J']:].flatten()[E]
                else:
                    b_stable[(len(species_params)**2 + (len(species_params)*k + j)*E.shape[0]):(len(species_params)**2 + (len(species_params)*k + j + 1)*E.shape[0])] = 0
        
    print('Third block of constructing A.')
    # Make constraints for partials of R_{k}^{F} \nu_{k}(F_{k}) and R_{k}^{ex} \nu_{k}^{ex}
    # with respect to x_{j}.
    for k in range(len(species_params)):
        for j in range(len(species_params)):
            A_kj_row_idx_start = int(len(species_params)**2 + (len(species_params)**2)*E.shape[0] + np.sum([len(species_params)*(len(species_params[w]['metabolites_idx'])+Jl) for w in range(k)]) + j*(len(species_params[k]['metabolites_idx'])+Jl))

            block_kj = S[:, lumen_reactions_idx.flatten()]
            block_kj = block_kj[np.concatenate([species_params[k]['metabolites_idx'].flatten(), lumen_metabolites_idx.flatten()]).flatten(), :]
            A_stable[A_kj_row_idx_start:(A_kj_row_idx_start + len(species_params[k]['metabolites_idx']) + Jl), (len(species_params)**2 + (len(species_params)*k + j)*Jl):(len(species_params)**2 + (len(species_params)*k + j + 1)*Jl)] = sparse.csr_matrix(block_kj)
            
            if len(F_vals[k]) == 0:
                pass
            else:
                block_kj_internal = S[:, species_params[k]['reactions_idx'][F_vals[k]].flatten()]
                block_kj_internal = block_kj_internal[np.concatenate([species_params[k]['metabolites_idx'], lumen_metabolites_idx]).flatten(), :]
                A_stable[A_kj_row_idx_start:(A_kj_row_idx_start + len(species_params[k]['metabolites_idx']) + Jl),
                        int((len(species_params)**2 + (len(species_params)**2)*Jl + len(species_params)*np.sum([F_vals[w].shape[0] for w in range(k)]) + j*F_vals[k].shape[0])):int((len(species_params)**2 + (len(species_params)**2)*Jl + len(species_params)*np.sum([F_vals[w].shape[0] for w in range(k)]) + (j+1)*F_vals[k].shape[0]))] = sparse.csr_matrix(block_kj_internal)
    
    #print('System of equations set up')
    result = sp_linalg.lsqr(A_stable, b_stable, atol=1e-6, btol=1e-6)[0]
    #print('System of equations solved')

    # Row k of h_vals gives values of partials of h_{k}.
    h_vals = np.zeros((len(species_params), len(species_params)))
    for k in range(len(species_params)):
        h_vals[k,:] = result[len(species_params)*k:len(species_params)*(k+1)]
    
    # Now construct the matrix M used to check stability.
    # Row k of M gives values of M_{k,j} for j = 0,1,2,3.
    M = np.zeros((len(species_params),len(species_params)))
    for k in range(len(species_params)):
        M[k,:] = bm_vals[k] * h_vals[k,:]
    
    # Compute eigenvalues of M. If real part of every eigenvalue
    # of M is less than 0, then we have a stable steady state.
    eigs = linalg.eig(M)[0]
    print('Eigenvalues computed')

    if np.max(eigs) > 1e-4:
        stable = False
    else:
        stable = True
    
    return (stable, np.max(eigs))
    
    

        