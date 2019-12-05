# sampler loop
# This version records everything
#
# Harness :: SeedOverRight , InitGibbState :: () -> NewState , transitionProposal :: (State-> StateChange)
#     , ApplytransitionUpdate :: (CurrentState, Proposal) -> NewState
#     , shouldItTerminate :: ( IterationCount,CurrentState,CurrentProposal -> Yes or NO  )
#  -> (finalState, Number of Steps, Initial rng state,Initial GibbsState Value, stateHistory, transitionProposalHistory)

harness <- function(seed=NA # seed is fine to leave as NA
    ,InitGibbState=NA
    ,TransitionProposal=NA
    ,ApplyTransition=NA
    ,ShouldWeTerminate=NA ) {
suppressWarnings({

if(is.na(InitGibbState)
    || is.na(TransitionProposal)
    || is.na(ApplyTransition)
    || is.na(ShouldWeTerminate)){
  errorCondition(message=
    "you must use:
      InitGibbState
      TransitionProposal
      ApplyTransition
      ShouldWeTerminate
    in keyword argument format
    "
    )
}
initialRngSnapshot <- tryCatch({.Random.seed},error=function(e){rnorm(1); .Random.seed})

restoreRNG <- FALSE

if (!is.na(seed))
  { set.seed(seed)
    restoreRNG <- TRUE
  }else {seed = initialRngSnapshot}

}) # end warning suppression

CurrentStep <- 0
CurrentState <- InitGibbState()
StateRecord <- c(CurrentState)
CurrentProposal <- TransitionProposal(CurrentState)
ProposalRecord <- c(CurrentProposal)

while(!ShouldWeTerminate(CurrentStep,CurrentState,CurrentProposal)){
  CurrentStep <- CurrentStep + 1
  CurrentState <- ApplyTransition(CurrentState,CurrentProposal)
  StateRecord <- c(StateRecord,CurrentState)
  CurrentProposal <- TransitionProposal(CurrentState)
  ProposalRecord <- c(ProposalRecord,CurrentProposal)
}

# restore RNG if required
if(restoreRNG)
  {set.seed(initialRngSnapshot)}

list(FinalState=CurrentState,StepCount=CurrentStep,TheStateRecord=StateRecord,TheProposalRecord=ProposalRecord,InitialRNGState=seed)
}

## Trivial example

# Gauss1dimWalk <- function(){
#   harness(
#     ,InitGibbState=function(){rnorm(1)}
#     ,TransitionProposal=function(x){rnorm(1) }# doesn't care about current state
#        # for graphs, you might depend on the # of nbrs of current vertex/state to determine how you sample over choosing which neighbor you propose moving to
#     ,ApplyTransition=function(state,proposal){state + proposal }
#     ,ShouldWeTerminate=function(step,state,proposal){ (step > 10 )|| sum(abs(proposal)) > 100 || sqrt(sum(state^2) > 400)
#     }
#     )
# }
# 
# Gauss1dimWalk()
