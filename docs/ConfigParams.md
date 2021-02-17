# JSON Config Paramters

Click the arrows on the overarching parameters to reveal sub-parameters and explanations for each one. Each field has a recommended 'default' value as per the sample config template provided in MEGneto/templates/sample_config.json, but please noet that these values may not work your analysis. 

<!--CONTACT-->

<details>
<summary>Contact</summary>
Where: the sendEmail function is called at the end of each fcp_# step, and contact is passed as a parameter
<br>
Meaning: Email address to which to send pipeline’s progress updates contained in square brackets (default = ["firstname.lastname@sickkids.ca"]).
</details>

<!--EPOCHING-->

<details>
  <summary>Epoching</summary>
  Where: fcp_1_TaskEpoching.
  <br>
  Meaning: Epoch the data into trials.
  <br><br>
    <ul>
      <!--EPOCHING.PERIOD--> 
      <li>
        <details>
          <summary>Period</summary>
          Where: fcp_1_RestingStateEpoching, line 97.
          <br>
          Meaning: Indicates epoch length for epoching resting state data (default = 30).
        </details>
      </li>
      <!--EPOCHING.TOTALTIME--> 
      <li>
        <details>
          <summary>Total time</summary>
          Where: Nowhere.
          <br>
          Meaning: Relic from an older resting state epoch strategy, will be deleted from config template.
        </details>
      </li>
     <!--EPOCHING.HEADMOTION--> 
      <li>
        <details>
          <summary>Head motion</summary>
          Where: fcp_1_TaskEpoching.
          <br>
          Meaning: Overall, the field within head motion indicates specifications for initial handling of head motion.
          <br><br>
            <ul>
              <!--EPOCHING.HEADMOTION.THRESHOLDING-->
              <li>
              <details>
                <summary>Threshold</summary>
                Where: fcp_1_TaskEpoching in head motion correction.
                <br>
                Meaning: Threshold for which to reject trials with head motion. This is based on the coil distances from the chosen reference location in mm. For examples, 5 indicates a threshold of 5mm, means movement within a trial from -4mm to 4mm in one direction gives a range of 8mm and thus this trial would be flagged (default = 10).
              </details>
            </li></ul>
        </details>
 </details>
 
 <!--CLEANING OPTIONS-->

<details>
  <summary>Cleaning Options</summary>
  Where: fcp_1_TaskEpoching.
  <br>
  Meaning: Overall, the fields within cleaning options specify the handling of artifacts, ICA, noisy trials, and bad channels.
  <br><br>
    <ul>
     <!--CLEANING OPTIONS.ARTIFACT--> 
      <li>
        <details>
          <summary>Artifact</summary>
          Where: fcp_1_TaskEpoching in artifact detection/rejection (all types).
          <br>
          Meaning: The fields within artifact specify how we want to deal with various parts involved in detecting and rejecting different artifacts.
          <br><br>
            <ul>
            <!--CLEANING OPTIONS.ARTIFACT.DETECTION--> 
            <li>
              <details>
                <summary>Detection</summary>
                Where: fcp_1_TaskEpoching in artifact detection.
                <br>
                Meaning: 0 or 1 (0=no, 1=yes) to indicate if we want to detect artifacts (default = 1).
              </details>
            </li>
            <!--CLEANING OPTIONS.ARTIFACT.REJECTION--> 
            <li>
              <details>
                <summary>Rejection</summary>
                Where: Nowhere: instead, in fcp_1_TaskEpoching under “Artifact Rejection”, there is a field “cfg.artfctdef.reject” that is equal to “complete”.
                <br>
                Meaning: Indicates how much we want to reject trials with artifacts (e.g. “complete” removes the entire trial).
              </details>
            </li>
            <!--CLEANING OPTIONS.ARTIFACT.MUSCLE--> 
            <li>
              <details>
                <summary>Muscle</summary>
                Where: fcp_1_TaskEpoching in artifact detection, setting up the cfg to call “ft_artifact_muscle”.
                <br>
                Meaning: Overall, the field within muscle specify how we want to deal with muscle artifacts.
              <br><br>
                <ul>
                  <!--CLEANING OPTIONS.ARTIFACT.MUSCLE.BPFILTER--> 
                  <li>
                    <details>
                      <summary>Bpfilter</summary>
                      Where: fcp_1_TaskEpoching in artifact detection, setting up the cfg to call “ft_artifact_muscle”.
                      <br>
                      Meaning: “Yes” or “no” to indicate whether or not we want to bandpass filter (default = "yes").
                    </details>
                  </li>
                  <!--CLEANING OPTIONS.ARTIFACT.MUSCLE.BPFREQ--> 
                  <li>
                    <details>
                      <summary>Bpfreq</summary>
                      Where: fcp_1_TaskEpoching in artifact detection, setting up the cfg to call “ft_artifact_muscle”. 
                      <br>
                      Meaning: [x,y] (x being lower band, y being higher band) to specify what frequency band the filter should be (default = [110,140]).
                    </details>
                  </li>
                  <!--CLEANING OPTIONS.ARTIFACT.MUSCLE.BPFILTORD--> 
                  <li>
                    <details>
                      <summary>Bpfiltord</summary>
                      Where: fcp_1_TaskEpoching in artifact detection, setting up the cfg to call “ft_artifact_muscle”.
                      <br>
                      Meaning: Specifies the fiter order (default = 8).
                    </details>
                  </li>
                  <!--CLEANING OPTIONS.ARTIFACT.MUSCLE.BPFILTTYPE--> 
                  <li>
                    <details>
                      <summary>Bpfilttype</summary>
                      Where: fcp_1_TaskEpoching in artifact detection, setting up the cfg to call “ft_artifact_muscle".
                      <br>
                      Meaning: Specifies the type of filter (default = “but” for butterworth).
                    </details>
                  </li>
                  <!--CLEANING OPTIONS.ARTIFACT.MUSCLE.HILBERT--> 
                  <li>
                    <details>
                      <summary>Hilbert</summary>
                      Where:fcp_1_TaskEpoching in artifact detection, setting up the cfg to call “ft_artifact_muscle".
                      <br>
                      Meaning: “Yes” or “no” to indicate if we want to perform a hilbert transform (default = "yes").
                    </details>
                  </li>
                  <!--CLEANING OPTIONS.ARTIFACT.MUSCLE.BOXCAR--> 
                  <li>
                    <details>
                      <summary>Boxcar</summary>
                      Where: fcp_1_TaskEpoching in artifact detection, setting up the cfg to call “ft_artifact_muscle”. 
                      <br>
                      Meaning: Specifies window length (time window) for the moving average filter. Also known as a boxcar car smoothing kernel or sliding average (default = 0.2).
                    </details>
                  </li>
                  <!--CLEANING OPTIONS.ARTIFACT.MUSCLE.CUTOFF--> 
                  <li>
                    <details>
                      <summary>Cutoff</summary>
                      Where: fcp_1_TaskEpoching in artifact detection, setting up the cfg to call “ft_artifact_muscle”.
                      <br>
                      Meaning: Specifies frequency at which to cut off the signal (default = 30). This depends on if you are interested in alpha or gamma band activity or a wider range as in evoked activity.
                    </details>
                  </li>
                  <!--CLEANING OPTIONS.ARTIFACT.MUSCLE.TRLPADDING--> 
                  <li>
                    <details>
                      <summary>Trlpadding</summary>
                      Where: fcp_1_TaskEpoching in artifact detection, setting up the cfg to call “ft_artifact_muscle”. 
                      <br>
                      Meaning: Allows data to be padded on either side of the trial with a specified length so that artifact detection/rejection are performed on those data segments (i.e. If you wish to include data prior to/post the trial are included). This field is measured in seconds (default = 0.5, which means you wish to pad with 0.5 seconds on the right and left side).  
                    </details>
                  </li>
                  <!--CLEANING OPTIONS.ARTIFACT.MUSCLE.FLTPADDING--> 
                  <li>
                    <details>
                      <summary>Fltpadding</summary>
                      Where: fcp_1_TaskEpoching in artifact detection, setting up the cfg to call “ft_artifact_muscle”. 
                      <br>
                      Meaning: Filter padding is needed because filters may cause edge effects detected in artifact-detection & mistaken for actual artifacts (only used for filtering, not artifact detection). This field is measured in seconds (default = 0.1, which means you wish to pad with 0.1 seconds on the right and left side). 
                    </details>
                  </li>
                  <!--CLEANING OPTIONS.ARTIFACT.MUSCLE.ARTPADDING--> 
                  <li>
                    <details>
                      <summary>Artpadding</summary>
                      Where: fcp_1_TaskEpoching in artifact detection, setting up the cfg to call “ft_artifact_muscle”. 
                      <br>
                      Meaning: Often, artifacts start/end a bit later than what is detected by the artifact detection system, thus artifact padding is used to extend the artifact timeperiod on either side. This field is measured in seconds (default = 0.1, which means you wish to pad with 0.1 seconds on the right and left side). 
                    </details>
                  </li>
                </ul>
              </details>
            </li> 
            <!--CLEANING OPTIONS.ARTIFACT.JUMP--> 
            <li>
              <details>
                <summary>Jump</summary>
                Where: fcp_1_TaskEpoching in artifact detection.
                <br>
                Meaning: Overall, the fields within jump specify how we want to deal with jump artifacts.
                <br><br>
                  <ul>
                    <!--CLEANING OPTIONS.ARTIFACT.JUMP.CUTOFF--> 
                    <li>
                      <details>
                        <summary>Cutoff</summary>
                        Where: fcp_1_TaskEpoching in artifact detection.
                        <br>
                        Meaning: Cutoff frequency indicating at what point the signal should be classified as a jump artifact (default = 35). This depends on if you are interested in alpha or gamma band activity or a wider range as in evoked activity. 
                      </details>
                    </li>
                </ul>
              </details>
            </li>
            <!--CLEANING OPTIONS.ARTIFACT.ICACLEAN--> 
            <li>
              <details>
                <summary>icaClean</summary>
                Where: fcp_2_PreprocessingICA before we do/don’t run ICA.
                <br>
                Meaning: 0 or 1 (0=no, 1=yes) to indicacte whether or not we want to perform ICA (default = 1).
              </details>
            </li>
            <!--CLEANING OPTIONS.ARTIFACT.rmNOISYTRIALS--> 
            <li>
              <details>
                <summary>rmNoisyTrials</summary>
                Where: fcp_2_PreprocessingICA. 
                <br>
                Meaning: 0 or 1 (0=no, 1=yes) to specify whether or not we want to remove noisy trials/artifacts (default = 1).
              </details>
            </li>
            <!--CLEANING OPTIONS.ARTIFACT.rmBADCHANNELS--> 
            <li>
              <details>
                <summary>rmBadChannels</summary>
                Where: fcp_3_ChannelRepair when checking if we want to remove channels.
                <br>
                Meaning: 0 or 1 (0=no, 1=yes) to indicate whether or not we want to remove bad channels (default = 1).
              </details>
            </li>
          </ul>
    </details>
 </details>
  
 <!--FILTERING PARAMETERS-->

<details>
  <summary>filteringParameters</summary>
  Where: fcp_2_PreprocessingICA in setting up cfg for ft_preprocessing.
  <br>
  Meaning: Overall, the fields within filteringParameters provide filtering specifications.
  <br><br>
    <ul>
      <!--FILTERING PARAMETERS.CHANNEL--> 
      <li>
        <details>
          <summary>Channel</summary>
          Where: fcp_2_PreprocessingICA in setting up cfg for ft_preprocessing.
          <br>
          Meaning: Specifies which data channels to look at: (1. MEG- replaced by all MEG channels, 2. MEGREF-replaced by all MEG reference channels, 3. REFGRAD, 4. REFMAG). default = [“MEG”, “MEGREF”, “REFGRAD”, “REFMAG”].
        </details>
      </li>
      <!--FILTERING PARAMETERS.DFTFILTER--> 
      <li>
        <details>
          <summary>Dftfilter</summary>
          Where: fcp_2_PreprocessingICA in setting up cfg for ft_preprocessing.
          <br>
          Meaning: “Yes” or “no” to indicate whether or not we want to apply a notch filter to the data to remove the 50Hz or 60Hz line noise components ('zeroing'). Default = "yes".
        </details>
      </li>
     <!--FILTERING PARAMETERS.DFTFREQ--> 
      <li>
        <details>
          <summary>Dftfreq</summary>
          Where: fcp_2_PreprocessingICA in setting up cfg for ft_preprocessing
          <br>
          Meaning: Indicates whether the power-line frequency to filter out is 50 or 60Hz and its harmonic frequency, which is next multiple of the power-line frequency (e.g. [60,120] or [50,100]).
        </details>
      </li>
      <!--FILTERING PARAMETERS.BPFILTER--> 
      <li>
        <details>
          <summary>Bpfilter</summary>
          Where: fcp_2_PreprocessingICA in setting up cfg for ft_preprocessing.
          <br>
          Meaning: “Yes” or “no” to indicate if we want to do a bandpass filter (default = "yes").
        </details>
      </li>
      <!--FILTERING PARAMETERS.BPFREQ--> 
      <li>
        <details>
          <summary>Bpfreq</summary>
          Where: fcp_2_PreprocessingICA in setting up cfg for ft_preprocessing.
          <br>
          Meaning: [x,y] to specify what frequency band the filter should be where x is the lower frequency band, and y is the higher frequency band (default = [1,150]).
        </details>
      </li>
      <!--FILTERING PARAMETERS.BPFILTORD--> 
      <li>
        <details>
          <summary>Bpfiltord</summary>
          Where: fcp_2_PreprocessingICA in setting up cfg for ft_preprocessing.
          <br>
          Meaning: Specifies the filter order (default = 8).
        </details>
      </li>
      <!--FILTERING PARAMETERS.SAMPLERATE--> 
      <li>
        <details>
          <summary>sampleRate</summary>
          Where: fcp_2_PreprocessingICA for downsampling data AND fcp_4_beamforming to resample the data.
          <br>
          Meaning: Rate at which data is sampled (how many data points per second, default = 300) .
        </details>
      </li>
      <!--FILTERING PARAMETERS.CTFLAYOUR--> 
      <li>
        <details>
          <summary>CTFlayout</summary>
          Where: End of fcp_2_5_checkpoint for displaying ica channels function.
          <br>
          Meaning: Indicates which MEG model you’re using (here, the CTF 151 model) so that it can plot results on a 2D image of the head with proper electrode positions (default = “CFT151.lay”). 
        </details>
      </li>
 </details>
        
 <!--TASKFUNC-->

<details>
  <summary>taskFunc</summary>
  Where: fcp_1_Task_Epoching for setting up the cfg for epoching.
  Meaning: Overall, the fields within taskFunc specify the parameters for epoching.
  <br><br>
    <ul>
      <!--TASKFUNC.FUNCTION-->
      <li>
        <details>
          <summary>Function</summary>
          Where: fcp_1_Task_Epoching for setting up the cfg for epoching.
          <br>
          Meaning: Name of a custom task epoching function to parse data into trials (designed for marker epoching). This will likely be “@searchTaskTrialFun” for you.
        </details>
      </li>
      <!--TASKFUNC.TYPE-->
      <li>
        <details>
          <summary>Type</summary>
          Where: fcp_1_Task_Epoching for setting up the cfg for epoching.
          Meaning: Indicates the type of function specified in the "Function" field. If using "@searchTaskTrialFun", this field should be entered as "anonymous" to indicate that it is an anonymous type of function (Matlab lingo, see details here: https://www.mathworks.com/help/matlab/matlab_prog/anonymous-functions.html).
        </details>
      </li>
     <!--TASKFUNC.FILE-->
      <li>
        <details>
          <summary>File</summary>
          Where: fcp_1_Task_Epoching for setting up the cfg for epoching.
        </details>
      </li>
     <!--TASKFUNC.WORKSPACE-->
      <li>
        <details>
          <summary>Workspace</summary>
          Where: fcp_1_Task_Epoching for setting up the cfg for epoching.
        </details>
      </li>
      <!--TASKFUNC.WITHINFILEPATH-->
      <li>
        <details>
          <summary>Within_file_path</summary>
          Where: fcp_1_Task_Epoching for setting up the cfg for epoching.
        </details>
      </li>
 </details>
  
 <!--TASK-->

<details>
  <summary>task</summary>
  Where: fcp_1_TaskEpoching and fcp_2_PreprocessingICA.
  Meaning: Overall, the fields within task help specify parameters for epoching and preprocessing.
  <br><br>
    <ul>
      <!--TASK.ISREST-->
      <li>
        <details>
          <summary>isRest</summary>
          Where: fcp_2_PreprocessingICA when we load subject specific data. 
          <br>
          Meaning: “0” or “1” to indicate whether or not we are dealing with rest data (default = 0). 
        </details>
      </li>
     <!--TASK.TRIALDEF--> 
      <li>
        <details>
          <summary>Trialdef</summary>
          Where: fcp_1_Task_Epoching for setting up the cfg for epoching (used in search TaskTrialFun in detail) and Fcp_5_freqanalysis and Fcp_5_task_Connectivity.
          <br>
          Meaning: Not assigned a specific value, rather it is the umbrella for the fields within it. Sets up parameters for trial definition. 
          <br><br>
            <ul>
              <!--TASK.TRIALDEF.DETAILS--> 
              <li>
              <details>
                <summary>Details</summary>
                Where: searchTaskTrialFun.m.
                <br>
                Meaning: Not assigned a specific value, rather it is the umbrella for the fields within it. 
                <br><br>
                <ul>
                  <!--TASK.TRIALDEF.DETAILS.NAME--> 
                  <li>
                  <details>
                    <summary>Name</summary>
                    Where: searchTaskTrialFun.m.
                    <br>  
                  </li>
                  <!--TASK.TRIALDEF.DETAILS.INCLUDEONCE--> 
                  <li>
                  <details>
                    <summary>includeOnce</summary>
                    Where: searchTaskTrialFun.m.
                    <br>
                    Meaning: Specifies additional events that should be included for stats generation. That is if you wish to generate reaction times between the stimulus presentation and markers included in this parameter. Currently must be left as [“”] as this functionality is not supported. 
                    <br>  
                  </li>
                  <!--TASK.TRIALDEF.DETAILS.EXCLUDE--> 
                  <li>
                  <details>
                    <summary>Exclude</summary>
                    Where: searchTaskTrialFun.m.
                    <br>
                    Meaning: Doesn't seem to be legitimately used. Leave as [“”]. 
                    <br>  
                  </li>
                  <!--TASK.TRIALDEF.DETAILS.INCLUDE--> 
                  <li>
                  <details>
                    <summary>Include</summary>
                    Where: searchTaskTrialFun.m.
                    <br>
                    Meaning: Specifies the marker you wish to gain reaction time information for (e.g. “Correct” means you will get the reaction time between stimulus presentation and a correct response).
                    <br>  
                  </li>
                  <!--TASK.TRIALDEF.DETAILS.COUNTONLY--> 
                  <li>
                  <details>
                    <summary>countOnly</summary>
                    Where: searchTaskTrialFun.m.
                    <br>
                    Meaning: “true” or “false” to specify whether you only want information on the number of trials (true for only want number of trials information, false for want more information including reaction times). Default = false.
                    <br>  
                  </li>
                </ul>
              </details>
            <!--TASK.TRIALDEF.LIGHT--> 
            <li>
              <details>
               <summary>Light</summary>
               Where: Nowhere.
               <br>
               Meaning: This is a relic from old code. Do not fill this field, it will be removed in future config templates.
               <ul>
                 <!--TASK.TRIALDEF.LIGHT.AVGSTARTTHRESH-->
                 <li>
                  <details>
                   <summary>avgstartThresh</summary>
                   Where: Nowhere.
                  </details>
                 </li>
                </ul>
              </details>
            <!--TASK.TRIALDEF.PARAMETERS--> 
            <li>
              <details>
               <summary>Parameters</summary>
               Where: searchTaskTrialFun and fcp_5_freqanalysis and fcp_5_taskconnectivity.
               <br>
               Meaning: Not assigned a specific value, rather it is the umbrella for the fields within it.
               <ul>
                 <!--TASK.TRIALDEF.PARAMETERS.T0SHIFT--> 
                 <li>
                  <details>
                   <summary>T0shift</summary>
                   Where: searchTaskTrialFun. 
                   <br>
                   Meaning: Specifies the delay time that needs to be corrected for in seconds (e.g. 0.023).
                  </details>
                 </li>
                 <!--TASK.TRIALDEF.PARAMETERS.TEPOCH--> 
                 <li>
                  <details>
                   <summary>tEpoch</summary>
                   Where: searchTaskTrialFun AND Fcp_5_freqanalysis AND Fcp_5_task_Connectivity when reshaping catmatrix into acceptable format.
                   <br>
                   Meaning: Specifies the epoch time window (default = [-2, 2]).
                  </details>
                 </li>
                </ul>
              </details>
             <!--TASK.TRIALDEF.MARKERS--> 
             <li>
              <details>
               <summary>Markers</summary>
               Where: fcp_1_TaskEpoching
               <br>
               Meaning: Not assigned a specific value, rather it is the umbrella for the fields within it.
               <!--TASK.TRIALDEF.MARKERS.CORRECT--> 
               <ul>
                 <li>
                  <details>
                   <summary>Correct</summary>
                   Where: searchTaskTrialFun.
                   <br>
                   Meaning: Specifies the markers that mark a correct trials (e.g. [“LeftCorrect”], [“RightCorrect”] or multiple such as [“LeftCorrect”, “RightCorrect”]).
                  </details>
                 </li>
                 <!--TASK.TRIALDEF.MARKERS.INCORRECT--> 
                 <li>
                  <details>
                   <summary>Incorrect</summary>
                   Where: Nowhere.
                   <br>
                   Meaning: This is a relic from old code. Do not fill in this field, it will be deleted from future config templates.
                  </details>
                 </li>
                 <!--TASK.TRIALDEF.MARKERS.T0MARKER--> 
                 <li>
                  <details>
                   <summary>T0marker</summary>
                   Where: fcp_1_TaskEpoching when t0 markers are grabbed (the plotTriggers function) AND search TaskTrialFun.
                   <br>
                   Meaning: Specifies the event type to epoch around. This will be the marker that define t=0 of each trigger, that is the presentation trigger (eg. “OfflneLightOn”, “LeftButtonPress”, etc.).
                  </details>
                 </li>
                 <!--TASK.TRIALDEF.MARKERS.NEWTRIG--> 
                 <li>
                  <details>
                   <summary>newTrig</summary>
                   Where: Nowhere
                   <br>
                   Meaning: This is a relic from old code. Do not fill in this field, it will be deleted from future config templates.
                  </details>
                 </li></ul>
              </details>
            </li></ul>
        </details>
 </details>
  
<!--BEAMFORMING-->

<details>
  <summary>Beamforming</summary>
  Where: Fcp_4_beamforming. 
  <br>
  Meaning: Not assigned a specific value, rather it is the umbrella for the fields within it. 
  <br><br>
    <ul>
      <!--BEAMFORMING.HEADMODEL--> 
      <li>
        <details>
          <summary>Headmodel</summary>
          Where: Fcp_4_beamforming. 
          <br>
          Meaning: Not assigned a specific value, rather it is the umbrella for the fields within it. 
          <ul>
          <!--BEAMFORMING.HEADMODEL.METHOD--> 
          <li>
            <details>
              <summary>Method</summary>
              Where: Fcp_4_beamforming when setting up cf to prepare the T1 head model AND participant specific head models. 
              <br>
              Meaning: Specifies what form the head model should be (default = "singleshell").
            </details>
          </li>
          <!--BEAMFORMING.HEADMODEL.UNITS--> 
          <li>
            <details>
              <summary>Units</summary>
              Where: Fcp_4_beamforming when setting up cf to prepare the T1 headmodel.
              <br>
              Meaning: Specifies units for the head model (default = "cm").
            </details>
         </li></ul>
        </details>
      </li>
       <!--BEAMFORMING.TEMPLATE--> 
      <li>
        <details>
          <summary>Template</summary>
          Where: Fcp_4_beamforming when constructing the grid for the T1 template model. 
          <br>
          Meaning: Not assigned a specific value, rather it is the umbrella for the fields within it. 
          <ul>
          <!--BEAMFORMING.TEMPLATE.GRID--> 
          <li>
            <details>
              <summary>Grid</summary>
              Where: Fcp_4_beamforming when constructing the grid for the T1 template model .
              <br>
              Meaning: Not assigned a specific value, rather it is the umbrella for the fields within it. 
              <ul>
                <!--BEAMFORMING.TEMPLATE.GRID.RESOLUTION--> 
                <li>
                  <details>
                    <summary>Resolution</summary>
                    Where: Fcp_4_beamforming when constructing the grid for the T1 template model. 
                    <br>
                    Meaning: Number of the resolution of the template MNI grid, defined in mm (default = 1).
                  </details>
                </li>
                <!--BEAMFORMING.TEMPLATE.TIGHT--> 
                <li>
                  <details>
                    <summary>Tight</summary>
                    Where: Fcp_4_beamforming when constructing the grid for the T1 template model. 
                    <br>
                    Meaning: "yes" or "no" (default = "yes"). 
                  </details>
               </li>
               <!--BEAMFORMING.TEMPLATE.TINWARDSHIFT--> 
                <li>
                  <details>
                    <summary>Inwardshift</summary>
                    Where: Fcp_4_beamforming when constructing the grid for the T1 template model. 
                    <br>
                    Meaning: Number that defines how much the innermost surface should be moved inward to constrain sources to be considered inside the source compartment (default = 0).
                  </details>
                </li></ul>
            </details>
          </li>
          <!--BEAMFORMING.TEMPLATE.COORDSYS--> 
          <li>
            <details>
              <summary>Coordsys</summary>
              Where: Fcp_4_beamforming when loading T1 template.
              <br>
              Meaning: The coordinate system that is used (default = “spm”).
            </details>
         </li></ul>
        </details>
      </li>
      <!--BEAMFORMING.ATLAS--> 
      <li>
        <details>
          <summary>Atlas</summary>
          Where: Fcp_4_beamforming.
          <br>
          Meaning: Not assigned a specific value, rather it is the umbrella for the fields within it. 
          <ul>
          <!--BEAMFORMING.ATLAS.FILEPATH--> 
          <li>
            <details>
              <summary>Filepath</summary>
              Where: Fcp_4_beamforming just after we perform actual beamforming. 
              <br>
              Meaning: Specifies filepath to the atlas we wish to project on (default = “/template/atlas/aal/ROI_MNI_V4.nii”).
            </details>
          </li>
          <!--BEAMFORMING.ATLAS.INPUTCOORD--> 
          <li>
            <details>
              <summary>Inputcoord</summary>
              Where: Fcp_4_beamforming for visualization of the T1 segmented head model (to check for alignment with grid).
              <br>
              Meaning: Default = “mni”.
            </details>
          </li></ul>
        </details>
       </li>
       <!--BEAMFORMING.CHECKMRIVOLUMES--> 
       <li>
        <details>
          <summary>checkMRIvolumes</summary>
          Where: Fcp_4_beamforming for visualization of participant’s segmented head model (to check for alignment).
          <br>
          Meaning: Not assigned a specific value, rather it is the umbrella for the fields within it. 
          <ul>
          <!--BEAMFORMING.CHECKMRIVOLUMES.METHOD--> 
          <li>
            <details>
              <summary>Method</summary>
              Where: Fcp_4_beamforming for visualization of participant’s segmented head model (to check for alignment). Only used if code for the visualization of segmented head model is uncommented.
              <br>
              Meaning: Specifies plotting method. Options are presented below though we tend to use “slice”: (1) “slice” (plots the data on a number of slices in the same plane). (2) “ortho” (plots the data on three orthogonal slices). (3) “surface” (plots the data on a 3D brain surface). (4) “glassbrain” (plots a max-projection through the brain). (5) “vertex” (plots the grid points or vertices scaled according to the functional value). (6) “cloud” (plot the data as clouds, spheres, or points scaled according to the functional value).
            </details>
          </li>
          <!--BEAMFORMING.CHECKMRIVOLUMES.SLIDESDIM--> 
          <li>
            <details>
              <summary>Slicesdim</summary>
              Where: Fcp_4_beamforming for visualization of participant’s segmented head model (to check for alignment). Only used if code for the visualization of segmented head model is uncommented.
              <br>
              Meaning: Only used when “method” is set to “slice”. This specifies the dimension to slice on. Options are: 1 (x-axis). 2 (y-axis). 3 (z-axis). Default = 3.
            </details>
          </li>
          <!--BEAMFORMING.CHECKMRIVOLUMES.NSLICES--> 
          <li>
            <details>
              <summary>Nslices</summary>
              Where: Fcp_4_beamforming for visualization of participant’s segmented head model (to check for alignment). Only used if code for the visualization of segmented head model is uncommented.
              <br>
              Meaning: Only used when “method” is set to “slice”. This will specify the number of slices (default = 20).
            </details>
          </li>
          <!--BEAMFORMING.CHECKMRIVOLUMES.FUNPARAMETER--> 
          <li>
            <details>
              <summary>Funparameter</summary>
              Where: Nowhere.
              <br>
              Meaning: Relic from old code. Do not fill this in, it will be deleted from future config templates.
          </li></ul>
        </details>
       </li>
       <!--BEAMFORMING.SUBJ--> 
       <li>
        <details>
          <summary>Subj</summary>
          Where: Fcp_4_beamforming to prepare the subject specific source model (with reference to the T1 template model).
          <br>
          Meaning: Not assigned a specific value, rather it is the umbrella for the fields within it. 
          <ul>
          <!--BEAMFORMING.SUBJ.GRID--> 
          <li>
            <details>
              <summary>Grid</summary>
              Where: Fcp_4_beamforming to prepare the subject specific source model (with reference to the T1 template model).
              <br>
              Meaning: Not assigned a specific value, rather it is the umbrella for the fields within it. 
              <ul>
              <!--BEAMFORMING.SUBJ.GRID.WARPMNI-->
              <li>
                <details>
                  <summary>Warpmni</summary>
                  Where: Fcp_4_beamforming to prepare the subject specific source model (with reference to the T1 template model).
                  <br>
                  Meaning: “Yes” or “no” to specify whether we want to warp the model to the T1 model (which acts as a control to normalize across all participants). Default = "yes".
                </details>
              </li>
              <!--BEAMFORMING.SUBJ.GRID.NONLINEAR-->
              <li>
                <details>
                  <summary>Nonlinear</summary>
                  Where: Fcp_4_beamforming to prepare the subject specific source model (with reference to the T1 template model).
                  <br>
                  Meaning: “Yes” or “no” to indicate whether non-linear normalization should be used (default = "yes").
                </details>
              </li>
              <!--BEAMFORMING.SUBJ.GRID.UNIT-->
              <li>
                <details>
                  <summary>Unit</summary>
                  Where: Fcp_4_beamforming to prepare the subject specific source model (with reference to the T1 template model). 
                  <br>
                  Meaning: Specify units (default = "cm").
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
          Where: Fcp_4_beamforming to compute leadfield. 
          <br>
          Meaning: Not assigned a specific value, rather it is the umbrella for the fields within it which specify parameters for computing the leadfield used for each participant. 
          <ul>
          <!--BEAMFORMING.LEADFIELD.REDUCERANK--> 
          <li>
            <details>
              <summary>Reducerank</summary>
              Where: Nowhere, but in “compute leadfield” there is a field “cfg.reducerank” that has an associated number.
              <br>
              Meaning: “Yes” or “no” to specify whether we want to reduce rank, which addresses depth bias in the computation of the forward model.
            </details>
           </li>
           <!--BEAMFORMING.LEADFIELD.NORMALIZE--> 
           <li>
            <details>
              <summary>Normalize</summary>
              Where: Fcp_4_beamforming to compute leadfield. 
              <br>
              Meaning: “Yes” or “no” to specify whether we want to normalize, which addresses depth bias in the computation of the forward model (default = "no").
            </details>
           </li></ul>
        </details>
       </li>
       <!--BEAMFORMING.TIMEDOMAIN--> 
       <li>
        <details>
          <summary>TimeDomain</summary>
          Where: Fcp_4_beamforming to compute the covariance matrix.
          <br>
          Meaning: Not assigned a specific value, rather it is the umbrella for the fields within it. 
          <ul>
           <!--BEAMFORMING.TIMEDOMAIN.COVARIANCE--> 
            <li>
              <details>
                <summary>Covariance</summary>
                Where: Fcp_4_beamforming to compute the covariance matrix.
                <br>
                Meaning: "Yes” or “no” to specify if a covariance matrix should be computed (default = "yes").
              </details>
             </li>
             <!--BEAMFORMING.TIMEDOMAIN.COVARIANCEWINDOW--> 
            <li>
              <details>
                <summary>Covariancewindow</summary>
                Where: Fcp_4_beamforming to compute the covariance matrix.
                <br>
                Meaning: Specifies window length for covariance matrix computation. Options include: [start end] in seconds or “all”’, “miniperiod”, “maxperiod”, “prestim”, “postim”. See documentation on ft_timelockanalysis (https://www.fieldtriptoolbox.org/reference/ft_timelockanalysis/) for detail. Default = "all".
              </details>
             </li>
             <!--BEAMFORMING.TIMEDOMAIN.VARTRLLENGTH--> 
            <li>
              <details>
                <summary>Vartrllength</summary>
                Where: fcp_4_beamforming prior to computing covariance matrix.
              </details>
             </li>
             <!--BEAMFORMING.TIMEDOMAIN.PROJECTMOM--> 
             <li>
              <details>
                <summary>Projectmom</summary>
                Where: Fcp_4_beamforming just after we perform actual beamforming.
                <br>
                Meaning: Projects data to the dominant eigenvector so we have one resulting timeseries per source reflecting overall change in activity. Else, there would be reconstructed source data in three orientations. Options are “yes” or “no” (default = "yes").
              </details>
             </li></ul>
        </details>
       </li>
       <!--BEAMFORMING.OPTIONS--> 
       <li>
        <details>
          <summary>Options</summary>
          Where: Fcp_4_beamforming.
          <br>
          Meaning: Not assigned a specific value, rather it is the umbrella for the fields within it. 
          <ul>
           <!--BEAMFORMING.OPTIONS.KEEPTRIALS--> 
            <li>
              <details>
                <summary>Keeptrials</summary>
                Where: Fcp_4_beamforming when computing the covariance matrix.
                <br>
                Meaning: Maintains separation of trials rather than returning the average source time series across all trials (default = "yes").
              </details>
             </li>
             <!--BEAMFORMING.OPTIONS.KEEPFILTER--> 
            <li>
              <details>
                <summary>Keepfilter</summary>
                Where: Fcp_4_beamforming after computing sensor weights.
                <br>
                Meaning: “Yes” or “no” to specify whether or not we want to keep the filer → which we do so that we can project all the data points through it to do the beamforming (default = "yes").
              </details>
             </li>
             <!--BEAMFORMING.OPTIONS.RAWTRIAL--> 
            <li>
              <details>
                <summary>Rawtrial</summary>
                Where: Fcp_4_beamforming just before we perform actual beamforming.
              </details>
            </li></ul>
        </details>
       </li>
       <!--BEAMFORMING.METHOD--> 
       <li>
        <details>
          <summary>Method</summary>
          Where: Fcp_4_beamforming just before we perform actual beamforming.
          <br>
          Meaning: Specifies what method we want to use to perform the beamforming (default = “lcmv”).
        </details>
      </li>
      <!--BEAMFORMING.REP_TIMESERIES-->
      <li>
        <details>
          <summary>Rep_timeseries</summary>
          Where: Fcp_4_beamforming just before a representative timeseries for each ROI. 
          <br>
          Meaning: Specifies what method we want to use to derive the timeseries. Options are "mean" or "pca" (default = "mean").
        </details>
      </li></ul>
 </details>
 
 <!--CONNECTIVITY-->

<details>
<summary>Connectivity</summary>
Where: Fcp_5_task_Connectivity.
Meaning: Not assigned a specific value, rather it is the umbrella for the fields within it. 
<br><br>
<ul>
  <!--CONNECTIVITY.METHOD-->
  <li>
    <details>
      <summary>Method</summary>
      Where: Fcp_5_task_Connectivity.
      <br>
      Meaning: Specifies what metric to use for connectivity analysis (default = "wpli_debiased").
    </details>
  </li>
  <!--CONNECTIVITY.FILTFREQS-->
  <li>
    <details>
      <summary>Filt_freqs</summary>
      Where: Fcp_5_task_Connectivity
      <br>
      Meaning: Specifies the various frequency bands (e.g. [4,7] for theta, [8,12] for alpha, etc. Default = [ [4,7], [8,12], [13,29],[30,59], [60,100] ].
    </details>
  </li>
  <!--CONNECTIVITY.COLLAPSEBAND-->
  <li>
    <details>
      <summary>Collapse_band</summary>
      Where: Nowhere.
      <br>
      Meaning: Specifies method of collapsing multiple reconstructed sources within a certain ROI to a representative timeseries for that ROI. Options are "mean" or "max" (default = "max").
    </details>
  </li></ul>
</details>

 
