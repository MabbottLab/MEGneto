# JSON Config Paramters

This file is the documentation for the paramters found in the JSON config file. Click the arrow on each parameter to reveal a dropdown explanation of where the paramter is used in the code and what it means. 

<!--CONTACT-->

<details>
<summary>Contact</summary>
Where: sendEmail function is called at the end of each fcp_# step, and contact is passed as a parameter.
<br>
Meaning: Email address to which to send pipeline’s progress updates (contained in square brackets, [ ])
</details>

<!--EPOCHING-->

<details>
  <summary>Epoching</summary>
  Where: fcp_1_TaskEpoching.
  <br>
  Meaning: Epoch the data into trials
  <br><br>
    <ul>
      <!--EPOCHING.PERIOD--> 
      <li>
        <details>
          <summary>Period</summary>
          Where: fcp_1_RestingStateEpoching, line 97
          <br>
          Meaning: Indicates epoch length for epoching resting state data
        </details>
      </li>
      <!--EPOCHING.TOTALTIME--> 
      <li>
        <details>
          <summary>Total time</summary>
          Where: Nowhere
          <br>
          Meaning: Relic from an older resting state epoch strategy. 
        </details>
      </li>
     <!--EPOCHING.HEADMOTION--> 
      <li>
        <details>
          <summary>Head motion</summary>
          Where: fcp_1_TaskEpoching
          <br>
          Meaning: Specifics for initial handling of head motion
          <br><br>
            <ul>
              <!--EPOCHING.HEADMOTION.THRESHOLDING-->
              <li>
              <details>
                <summary>Threshold</summary>
                Where: fcp_1_TaskEpoching in head motion correction
                <br>
                Meaning: Threshold for which to reject trials with head motion 
              </details>
            </li></ul>
        </details>
 </details>
 
 <!--CLEANING OPTIONS-->

<details>
  <summary>Cleaning Options</summary>
  Where: fcp_1_TaskEpoching
  <br>
  Meaning: Overall, gives options for how to handle artifacts, ICA, noisy trials, bad channels
  <br><br>
    <ul>
     <!--CLEANING OPTIONS.ARTIFACT--> 
      <li>
        <details>
          <summary>Artifact</summary>
          Where: fcp_1_TaskEpoching in artifact detection/rejection (all types)
          <br>
          Meaning: Specifies how we want to deal with various parts involved in detecting and rejecting different artifacts
          <br><br>
            <ul>
            <!--CLEANING OPTIONS.ARTIFACT.DETECTION--> 
            <li>
              <details>
                <summary>Detection</summary>
                Where: fcp_1_TaskEpoching in artifact detection
                <br>
                Meaning: “0” or “1” to indicate if we want to detect artifacts
              </details>
            </li>
            <!--CLEANING OPTIONS.ARTIFACT.REJECTION--> 
            <li>
              <details>
                <summary>Rejection</summary>
                Where: Nowhere: instead, in fcp_1_TaskEpoching under “Artifact Rejection”, there is a field “cfg.artfctdef.reject” that is equal to “complete” 
                <br>
                Meaning: Indicates how much we want to reject trials with artifacts (e.g. “complete” removes the entire trial)
              </details>
            </li>
            <!--CLEANING OPTIONS.ARTIFACT.MUSCLE--> 
            <li>
              <details>
                <summary>Muscle</summary>
                Where: fcp_1_TaskEpoching in artifact detection, setting up the cfg to call “ft_artifact_muscle”
                <br>
                Meaning: Configuration for muscle specifies how we want to deal with muscle artifacts 
              <br><br>
                <ul>
                  <!--CLEANING OPTIONS.ARTIFACT.MUSCLE.BPFILTER--> 
                  <li>
                    <details>
                      <summary>Bpfilter</summary>
                      Where: fcp_1_TaskEpoching in artifact detection, setting up the cfg to call “ft_artifact_muscle” 
                      <br>
                      Meaning: “Yes” or “no” to indicate whether or not we want to bandpass filter
                    </details>
                  </li>
                  <!--CLEANING OPTIONS.ARTIFACT.MUSCLE.BPFREQ--> 
                  <li>
                    <details>
                      <summary>Bpfreq</summary>
                      Where: fcp_1_TaskEpoching in artifact detection, setting up the cfg to call “ft_artifact_muscle” 
                      <br>
                      Meaning: [x,y] to specify what frequency band the filter should be
                    </details>
                  </li>
                  <!--CLEANING OPTIONS.ARTIFACT.MUSCLE.BPFILTORD--> 
                  <li>
                    <details>
                      <summary>Bpfiltord</summary>
                      Where: fcp_1_TaskEpoching in artifact detection, setting up the cfg to call “ft_artifact_muscle”
                      <br>
                      Meaning: Specifies the fiter orde
                    </details>
                  </li>
                  <!--CLEANING OPTIONS.ARTIFACT.MUSCLE.BPFILTTYPE--> 
                  <li>
                    <details>
                      <summary>Bpfilttype</summary>
                      Where: fcp_1_TaskEpoching in artifact detection, setting up the cfg to call “ft_artifact_muscle
                      <br>
                      Meaning: Specifies the type of filter (e.g. “but” for butterworth)
                    </details>
                  </li>
                  <!--CLEANING OPTIONS.ARTIFACT.MUSCLE.HILBERT--> 
                  <li>
                    <details>
                      <summary>Hilbert</summary>
                      Where:fcp_1_TaskEpoching in artifact detection, setting up the cfg to call “ft_artifact_muscle
                      <br>
                      Meaning: “Yes” or “no” to indicate if we want to perform a hilbert transform
                    </details>
                  </li>
                  <!--CLEANING OPTIONS.ARTIFACT.MUSCLE.BOXCAR--> 
                  <li>
                    <details>
                      <summary>Boxcar</summary>
                      Where: fcp_1_TaskEpoching in artifact detection, setting up the cfg to call “ft_artifact_muscle” 
                      <br>
                      Meaning: Specifies window length for the moving average filter. Also known as a boxcar car smoothing kernel or sliding average (aka window length
                    </details>
                  </li>
                  <!--CLEANING OPTIONS.ARTIFACT.MUSCLE.CUTOFF--> 
                  <li>
                    <details>
                      <summary>Cutoff</summary>
                      Where: fcp_1_TaskEpoching in artifact detection, setting up the cfg to call “ft_artifact_muscle”
                      <br>
                      Meaning: Specifies frequency at which to cut off the signal
                    </details>
                  </li>
                  <!--CLEANING OPTIONS.ARTIFACT.MUSCLE.TRLPADDING--> 
                  <li>
                    <details>
                      <summary>Trlpadding</summary>
                      Where: fcp_1_TaskEpoching in artifact detection, setting up the cfg to call “ft_artifact_muscle” 
                      <br>
                      Meaning: Allows data to be padded on either side of the trial with a specified length so that artifact detection/rejection are performed on those data segments (i.e. If you wish to include data prior to/post the trial are included) 
                    </details>
                  </li>
                  <!--CLEANING OPTIONS.ARTIFACT.MUSCLE.FLTPADDING--> 
                  <li>
                    <details>
                      <summary>Fltpadding</summary>
                      Where: fcp_1_TaskEpoching in artifact detection, setting up the cfg to call “ft_artifact_muscle” 
                      <br>
                      Meaning: Only used for filtering, not artifact detection. Since filters may cause edge effects detected in artifact-detection & mistaken for actual artifacts, we need filter padding. This reads additional data on either side before filtering.
                    </details>
                  </li>
                  <!--CLEANING OPTIONS.ARTIFACT.MUSCLE.ARTPADDING--> 
                  <li>
                    <details>
                      <summary>Artpadding</summary>
                      Where: fcp_1_TaskEpoching in artifact detection, setting up the cfg to call “ft_artifact_muscle” 
                      <br>
                      Meaning: Often, artifacts start/end a bit later than what is detected by the artifact detection system. Thus artifact padding is used to extend the artifact timeperiod on either side.
                    </details>
                  </li>
                </ul>
              </details>
            </li> 
            <!--CLEANING OPTIONS.ARTIFACT.JUMP--> 
            <li>
              <details>
                <summary>Jump</summary>
                Where: fcp_1_TaskEpoching in artifact detection
                <br>
                Meaning: Specifies how we want to deal with jump artifacts
                <br><br>
                  <ul>
                    <!--CLEANING OPTIONS.ARTIFACT.JUMP.CUTOFF--> 
                    <li>
                      <details>
                        <summary>Cutoff</summary>
                        Where: fcp_1_TaskEpoching in artifact detection
                        <br>
                        Meaning: Cutoff frequency indicating at what point the signal should be classified as a jump artifact 
                      </details>
                    </li>
                </ul>
              </details>
            </li>
            <!--CLEANING OPTIONS.ARTIFACT.ICACLEAN--> 
            <li>
              <details>
                <summary>icaClean</summary>
                Where: fcp_2_PreprocessingICA before we do/don’t run ICA
                <br>
                Meaning: “0” or “1” to indicacte whether or not we want to perform ICA
              </details>
            </li>
            <!--CLEANING OPTIONS.ARTIFACT.rmNOISYTRIALS--> 
            <li>
              <details>
                <summary>rmNoisyTrials</summary>
                Where: fcp_2_PreprocessingICA 
                <br>
                Meaning: “0” or “1” to specify whether or not we want to remove noisy trials (artifacts)
              </details>
            </li>
            <!--CLEANING OPTIONS.ARTIFACT.rmBADCHANNELS--> 
            <li>
              <details>
                <summary>rmBadChannels</summary>
                Where: fcp_3_ChannelRepair when checking if we want to remove channels
                <br>
                Meaning: “0” or “1” to indicate whether or not we want to remove bad channels (which we kept a list of from fcp_1)
              </details>
            </li>
          </ul>
    </details>
 </details>
  
 <!--FILTERING PARAMETERS-->

<details>
  <summary>filteringParameters</summary>
  Where: fcp_2_PreprocessingICA in setting up cfg for ft_preprocessing
  <br>
  Meaning: Overall, provides filtering specifications
  <br><br>
    <ul>
      <!--FILTERING PARAMETERS.CHANNEL--> 
      <li>
        <details>
          <summary>Channels</summary>
          Where: fcp_2_PreprocessingICA in setting up cfg for ft_preprocessing
          <br>
          Meaning: Specifies which data channels to look at: (1. MEG- replaced by all MEG channels, 2. MEGREF-replaced by all MEG reference channels, 3. REFGRAD, 4. REFMAG)
        </details>
      </li>
      <!--FILTERING PARAMETERS.DFTFILTER--> 
      <li>
        <details>
          <summary>Dftfilter</summary>
          Where: fcp_2_PreprocessingICA in setting up cfg for ft_preprocessing
          <br>
          Meaning: “Yes” or “no” to indicate whether or not we want to apply a notch filter to the data to remove the 50Hz
 or 60Hz line noise components ('zeroing').
        </details>
      </li>
     <!--FILTERING PARAMETERS.DFTFREQ--> 
      <li>
        <details>
          <summary>Dftfreq</summary>
          Where: fcp_2_PreprocessingICA in setting up cfg for ft_preprocessing
          <br>
          Meaning: Indicates whether the frequency to filter out is 50 or 60Hz. 
        </details>
      </li>
      <!--FILTERING PARAMETERS.BPFILTER--> 
      <li>
        <details>
          <summary>Bpfilter</summary>
          Where: fcp_2_PreprocessingICA in setting up cfg for ft_preprocessing
          <br>
          Meaning: “Yes” or “not” to indicate if we want to do a bandpass filter 
        </details>
      </li>
      <!--FILTERING PARAMETERS.BPFREQ--> 
      <li>
        <details>
          <summary>Bpfreq</summary>
          Where: fcp_2_PreprocessingICA in setting up cfg for ft_preprocessing
          <br>
          Meaning: [x,y] to specify what frequency band the filter should be
        </details>
      </li>
      <!--FILTERING PARAMETERS.BPFILTORD--> 
      <li>
        <details>
          <summary>Bpfiltord</summary>
          Where: fcp_2_PreprocessingICA in setting up cfg for ft_preprocessing
          <br>
          Meaning: Specifies the filter order
        </details>
      </li>
      <!--FILTERING PARAMETERS.SAMPLERATE--> 
      <li>
        <details>
          <summary>sampleRate</summary>
          Where: fcp_2_PreprocessingICA for downsampling data AND fcp_4_beamforming to resample the data
          <br>
          Meaning: Rate at which data is sampled (how many data points per second) 
        </details>
      </li>
      <!--FILTERING PARAMETERS.CTFLAYOUR--> 
      <li>
        <details>
          <summary>CTFlayout</summary>
          Where: End of fcp_2_5_checkpoint for displaying ica channels functio
          <br>
          Meaning: Indicates which MEG model you’re using (here, the CTF 151 model) so that it can plot results on a 2D image of the head with proper electrode positions.
        </details>
      </li>
 </details>
        
 <!--TASKFUNC-->

<details>
  <summary>taskFunc</summary>
  Where: fcp_1_Task_Epoching for setting up the cfg for epoching
  <br><br>
    <ul>
      <!--TASKFUNC.FUNCTION-->
      <li>
        <details>
          <summary>Function</summary>
          Where: fcp_1_Task_Epoching for setting up the cfg for epoching
          <br>
          Meaning: Name of a custom task epoching function to parse data into trials. Designed for marker epoching.
        </details>
      </li>
      <!--TASKFUNC.TYPE-->
      <li>
        <details>
          <summary>Type</summary>
          Where: fcp_1_Task_Epoching for setting up the cfg for epoching
        </details>
      </li>
     <!--TASKFUNC.FILE-->
      <li>
        <details>
          <summary>File</summary>
          Where: fcp_1_Task_Epoching for setting up the cfg for epoching
        </details>
      </li>
     <!--TASKFUNC.WORKSPACE-->
      <li>
        <details>
          <summary>Workspace</summary>
          Where: fcp_1_Task_Epoching for setting up the cfg for epoching
        </details>
      </li>
      <!--TASKFUNC.WITHINFILEPATH-->
      <li>
        <details>
          <summary>Within_file_path</summary>
          Where: fcp_1_Task_Epoching for setting up the cfg for epoching
        </details>
      </li>
 </details>
  
 <!--TASK-->

<details>
  <summary>task</summary>
  Where: fcp_1_TaskEpoching and fcp_2_PreprocessingICA
  <br><br>
    <ul>
      <!--TASK.ISREST-->
      <li>
        <details>
          <summary>isRest</summary>
          Where: fcp_2_PreprocessingICA when we load subject specific data 
          <br>
          Meaning: “0” or “1” to indicate whether or not we are dealing with rest data 
        </details>
      </li>
     <!--TASK.TRIALDEF--> 
      <li>
        <details>
          <summary>Trialdef</summary>
          Where: fcp_1_Task_Epoching for setting up the cfg for epoching (used in search TaskTrialFun in detail) and Fcp_5_freqanalysis and Fcp_5_task_Connectivity
          <br><br>
            <ul>
              <!--TASK.TRIALDEF.DETAILS--> 
              <li>
              <details>
                <summary>Details</summary>
                Where: searchTaskTrialFun
                <br><br>
                <ul>
                  <!--TASK.TRIALDEF.DETAILS.NAME--> 
                  <li>
                  <details>
                    <summary>Name</summary>
                    Where: searchTaskTrialFun
                    <br>  
                  </li>
                  <!--TASK.TRIALDEF.DETAILS.INCLUDEONCE--> 
                  <li>
                  <details>
                    <summary>includeOnce</summary>
                    Where: searchTaskTrialFun
                    <br>  
                  </li>
                  <!--TASK.TRIALDEF.DETAILS.EXCLUDE--> 
                  <li>
                  <details>
                    <summary>Exclude</summary>
                    Where: searchTaskTrialFun
                    <br>  
                  </li>
                  <!--TASK.TRIALDEF.DETAILS.INCLUDE--> 
                  <li>
                  <details>
                    <summary>Include</summary>
                    Where: searchTaskTrialFun
                    <br>  
                  </li>
                  <!--TASK.TRIALDEF.DETAILS.COUNTONLY--> 
                  <li>
                  <details>
                    <summary>countOnly</summary>
                    Where: searchTaskTrialFun
                    <br>  
                  </li>
                </ul>
              </details>
            <!--TASK.TRIALDEF.LIGHT--> 
            <li>
              <details>
               <summary>Light</summary>
               Where: Nowhere
               <ul>
                 <!--TASK.TRIALDEF.LIGHT.AVGSTARTTHRESH-->
                 <li>
                  <details>
                   <summary>avgstartThresh</summary>
                   Where: Nowhere
                  </details>
                 </li>
                </ul>
              </details>
            <!--TASK.TRIALDEF.PARAMETERS--> 
            <li>
              <details>
               <summary>Parameters</summary>
               Where: searchTaskTrialFun and fcp_5_freqanalysis and fcp_5_taskconnectivity
               <ul>
                 <!--TASK.TRIALDEF.PARAMETERS.T0SHIFT--> 
                 <li>
                  <details>
                   <summary>T0shift</summary>
                   Where: searchTaskTrialFun 
                  </details>
                 </li>
                 <!--TASK.TRIALDEF.PARAMETERS.TEPOCH--> 
                 <li>
                  <details>
                   <summary>tEpoch</summary>
                   Where: searchTaskTrialFun AND Fcp_5_freqanalysis AND Fcp_5_task_Connectivity when reshaping catmatrix into acceptable format 
                  </details>
                 </li>
                </ul>
              </details>
             <!--TASK.TRIALDEF.MARKERS--> 
             <li>
              <details>
               <summary>Markers</summary>
               Where: fcp_1_TaskEpoching
               <!--TASK.TRIALDEF.MARKERS.CORRECT--> 
               <ul>
                 <li>
                  <details>
                   <summary>Correct</summary>
                   Where: searchTaskTrialFun 
                  </details>
                 </li>
                 <!--TASK.TRIALDEF.MARKERS.INCORRECT--> 
                 <li>
                  <details>
                   <summary>Incorrect</summary>
                   Where: Nowhere
                  </details>
                 </li>
                 <!--TASK.TRIALDEF.MARKERS.T0MARKER--> 
                 <li>
                  <details>
                   <summary>T0marker</summary>
                   Where: fcp_1_TaskEpoching when t0 markers are grabbed (the plotTriggers function) AND search TaskTrialFun
                   <br>
                   Meaning: Specifies the event type to epoch around (eg. OfflneLightOn, LeftButtonPress, etc.)
                  </details>
                 </li>
                 <!--TASK.TRIALDEF.MARKERS.NEWTRIG--> 
                 <li>
                  <details>
                   <summary>newTrig</summary>
                   Where: Nowhere
                  </details>
                 </li></ul>
              </details>
            </li></ul>
        </details>
 </details>
  
<!--BEAMFORMING-->

<details>
  <summary>Beamforming</summary>
  Where: Fcp_4_beamforming 
  <br><br>
    <ul>
      <!--BEAMFORMING.HEADMODEL--> 
      <li>
        <details>
          <summary>Headmodel</summary>
          Where: Fcp_4_beamforming 
          <ul>
          <!--BEAMFORMING.HEADMODEL.METHOD--> 
          <li>
            <details>
              <summary>Method</summary>
              Where: Fcp_4_beamforming when setting up cf to prepare the T1 head model AND participant specific head models 
              <br>
              Meaning: Specifies what form the head model should be (e.g. single shell)
            </details>
          </li>
          <!--BEAMFORMING.HEADMODEL.UNITS--> 
          <li>
            <details>
              <summary>Units</summary>
              Where: Fcp_4_beamforming when setting up cf to prepare the T1 headmodel
              <br>
              Meaning: Specifies units for the head model (e.g. cm)
            </details>
         </li></ul>
        </details>
      </li>
       <!--BEAMFORMING.TEMPLATE--> 
      <li>
        <details>
          <summary>Template</summary>
          Where: Fcp_4_beamforming when constructing the grid for the T1 template model 
          <ul>
          <!--BEAMFORMING.TEMPLATE.GRID--> 
          <li>
            <details>
              <summary>Grid</summary>
              Where: Fcp_4_beamforming when constructing the grid for the T1 template model 
              <ul>
                <!--BEAMFORMING.TEMPLATE.GRID.RESOLUTION--> 
                <li>
                  <details>
                    <summary>Resolution</summary>
                    Where: Fcp_4_beamforming when constructing the grid for the T1 template model 
                  </details>
                </li>
                <!--BEAMFORMING.TEMPLATE.TIGHT--> 
                <li>
                  <details>
                    <summary>Tight</summary>
                    Where: Fcp_4_beamforming when constructing the grid for the T1 template model 
                  </details>
               </li>
               <!--BEAMFORMING.TEMPLATE.TINWARDSHIFT--> 
                <li>
                  <details>
                    <summary>Inwardshift</summary>
                    Where: Fcp_4_beamforming when constructing the grid for the T1 template model 
                  </details>
                </li></ul>
            </details>
          </li>
          <!--BEAMFORMING.TEMPLATE.COORDSYS--> 
          <li>
            <details>
              <summary>Coordsys</summary>
              Where: Fcp_4_beamforming when loading T1 template
            </details>
         </li></ul>
        </details>
      </li>
      <!--BEAMFORMING.ATLAS--> 
      <li>
        <details>
          <summary>Atlas</summary>
          Where: Fcp_4_beamforming
          <ul>
          <!--BEAMFORMING.ATLAS.FILEPATH--> 
          <li>
            <details>
              <summary>Filepath</summary>
              Where: Fcp_4_beamforming just after we perform actual beamforming 
              <br>
              Meaning: Specifies filepath to the atlas we wish to project on
            </details>
          </li>
          <!--BEAMFORMING.ATLAS.INPUTCOORD--> 
          <li>
            <details>
              <summary>Inputcoord</summary>
              Where: Fcp_4_beamforming for visualization of the T1 segmented head model (to check for alignment with grid) 
            </details>
          </li></ul>
        </details>
       </li>
       <!--BEAMFORMING.CHECKMRIVOLUMES--> 
       <li>
        <details>
          <summary>checkMRIvolumes</summary>
          Where: Fcp_4_beamforming for visualization of participant’s segmented head model (to check for alignment)
          <ul>
          <!--BEAMFORMING.CHECKMRIVOLUMES.METHOD--> 
          <li>
            <details>
              <summary>Method</summary>
              Where: Fcp_4_beamforming for visualization of participant’s segmented head model (to check for alignment)
            </details>
          </li>
          <!--BEAMFORMING.CHECKMRIVOLUMES.SLIDESDIM--> 
          <li>
            <details>
              <summary>Slidesdim</summary>
              Where: Fcp_4_beamforming for visualization of participant’s segmented head model (to check for alignment)
            </details>
          </li>
          <!--BEAMFORMING.CHECKMRIVOLUMES.NSLICES--> 
          <li>
            <details>
              <summary>Nslices</summary>
              Where: Fcp_4_beamforming for visualization of participant’s segmented head model (to check for alignment)
            </details>
          </li>
          <!--BEAMFORMING.CHECKMRIVOLUMES.FUNPARAMETER--> 
          <li>
            <details>
              <summary>Funparameter</summary>
              Where: Nowhere
              <br>
              Meaning: Functional parameter set to be ‘power’, which indicates what to plot on the anatomical data.
            </details>
          </li></ul>
        </details>
       </li>
       <!--BEAMFORMING.SUBJ--> 
       <li>
        <details>
          <summary>Subj</summary>
          Where: Fcp_4_beamforming to prepare the subject specific source model (with reference to the T1 template model)
          <ul>
          <!--BEAMFORMING.SUBJ.GRID--> 
          <li>
            <details>
              <summary>Grid</summary>
              Where: Fcp_4_beamforming to prepare the subject specific source model (with reference to the T1 template model)
              <ul>
              <!--BEAMFORMING.SUBJ.GRID.WARPMNI-->
              <li>
                <details>
                  <summary>Warpmni</summary>
                  Where: Fcp_4_beamforming to prepare the subject specific source model (with reference to the T1 template model)
                  <br>
                  Meaning: “Yes” or “no” to specify whether we want to warp the model to the T1 model (which acts as a control to normalize across all participants) 
                </details>
              </li>
              <!--BEAMFORMING.SUBJ.GRID.NONLINEAR-->
              <li>
                <details>
                  <summary>Nonlinear</summary>
                  Where: Fcp_4_beamforming to prepare the subject specific source model (with reference to the T1 template model)
                </details>
              </li>
              <!--BEAMFORMING.SUBJ.GRID.UNIT-->
              <li>
                <details>
                  <summary>Unit</summary>
                  Where: Fcp_4_beamforming to prepare the subject specific source model (with reference to the T1 template model)
                  <br>
                  Meaning: Specify units (e.g. cm)
                </details>
              </li></ul>
            </details>
          </li></ul>
        </details>
       </li>
       <!--BEAMFORMING.LEADFIELD--> 
       <li>
        <details>
          <summary>Leadfield</summary>
          Where: Fcp_4_beamforming to compute leadfield 
          <br>
          Meaning: Computes leadfield used for each participant 
          <ul>
          <!--BEAMFORMING.LEADFIELD.REDUCERANK--> 
          <li>
            <details>
              <summary>Reducerank</summary>
              Where: Nowhere, BUT in “compute leadfield” there is a field “cfg.reducerank” that has an associated number
              <br>
              Meaning: “Yes” or “no” to specify whether we want to reduce rank. Addresses depth bias in the computation of the forward model.
            </details>
           </li>
           <!--BEAMFORMING.LEADFIELD.NORMALIZE--> 
           <li>
            <details>
              <summary>Normalize</summary>
              Where: Fcp_4_beamforming to compute leadfield 
              <br>
              Meaning: “Yes” or “no” to specify whether we want to normalize. Addresses depth bias in the computation of the forward model. 
            </details>
           </li></ul>
        </details>
       </li>
       <!--BEAMFORMING.TIMEDOMAIN--> 
       <li>
        <details>
          <summary>TimeDomain</summary>
          Where: Fcp_4_beamforming to compute the covariance matrix
          <ul>
           <!--BEAMFORMING.TIMEDOMAIN.COVARIANCE--> 
            <li>
              <details>
                <summary>Covariance</summary>
                Where: Fcp_4_beamforming to compute the covariance matrix
              </details>
             </li>
             <!--BEAMFORMING.TIMEDOMAIN.COVARIANCEWINDOW--> 
            <li>
              <details>
                <summary>Covariancewindow</summary>
                Where: Fcp_4_beamforming to compute the covariance matrix
              </details>
             </li>
             <!--BEAMFORMING.TIMEDOMAIN.VARTRLLENGTH--> 
            <li>
              <details>
                <summary>Vartrllength</summary>
                Where: Fcp_4_beamforming to compute the covariance matrix
              </details>
             </li>
             <!--BEAMFORMING.TIMEDOMAIN.PROJECTMOM--> 
             <li>
              <details>
                <summary>Projectmom</summary>
                Where: Fcp_4_beamforming just after we perform actual beamforming 
                <br>
                Meaning: Projects data to the dominant eigenvector so we have one resulting timeseries per source reflecting overall change in activity. Else, there would be reconstructed source data in three orientations.
              </details>
             </li></ul>
        </details>
       </li>
       <!--BEAMFORMING.OPTIONS--> 
       <li>
        <details>
          <summary>Options</summary>
          Where: Fcp_4_beamforming
          <ul>
           <!--BEAMFORMING.OPTIONS.KEEPTRIALS--> 
            <li>
              <details>
                <summary>Keeptrials</summary>
                Where: Fcp_4_beamforming when computing the covariance matrix
                <br>
                Meaning: Maintains separation of trials rather than returning the average source time series across all trials.
              </details>
             </li>
             <!--BEAMFORMING.OPTIONS.KEEPFILTER--> 
            <li>
              <details>
                <summary>Keepfilter</summary>
                Where: Fcp_4_beamforming after computing sensor weights
                <br>
                Meaning: “Yes” or “no” to specify whether or not we want to keep the filer → which we do so that we can project all the data points through it to do the beamforming
              </details>
             </li>
             <!--BEAMFORMING.OPTIONS.RAWTRIAL--> 
            <li>
              <details>
                <summary>Rawtrial</summary>
                Where: Fcp_4_beamforming just before we perform actual beamforming
              </details>
            </li></ul>
        </details>
       </li>
       <!--BEAMFORMING.METHOD--> 
       <li>
        <details>
          <summary>Method</summary>
          Where:Fcp_4_beamforming just before we perform actual beamforming 
          <br>
          Meaning: Specifies what method we want to use to perform the beam forming 
        </details>
      </li></ul>
 </details>
 
 <!--CONNECTIVITY-->

<details>
<summary>Connectivity</summary>
Where: Fcp_5_task_Connectivity
<br><br>
<ul>
  <!--CONNECTIVITY.METHOD-->
  <li>
    <details>
      <summary>Method</summary>
      Where: Fcp_5_task_Connectivity
      <br>
      Meaning: Specifies what metric to use for connectivity analysis (e.g. wpli, etc.)
    </details>
  </li>
  <!--CONNECTIVITY.FILTFREQS-->
  <li>
    <details>
      <summary>Filt_freqs</summary>
      Where: Fcp_5_task_Connectivity
      <br>
      Meaning: Specifies the various frequency bands (theta, alpha, beta, gamma, hi-gamma)
    </details>
  </li>
  <!--CONNECTIVITY.COLLAPSEBAND-->
  <li>
    <details>
      <summary>Collapse_band</summary>
      Where: Nowhere
      <br>
      Meaning: Specifies method of collapsing multiple reconstructed sources within a certain ROI to a representative timeseries for that ROI (e.g. mean or max)
    </details>
  </li></ul>
</details>

 
