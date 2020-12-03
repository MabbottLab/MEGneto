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
              </details>
            </li> 
            <!--CLEANING OPTIONS.ARTIFACT.JUMP--> 
            <li>
              <details>
                <summary>Jump</summary>
                Where: fcp_1_TaskEpoching in artifact detection
                <br>
                Meaning: Specifies how we want to deal with jump artifacts
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
