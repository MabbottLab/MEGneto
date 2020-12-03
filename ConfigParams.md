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
    
  <li>
    <details>
      <summary>Total time</summary>
      Where: Nowhere
      <br>
      Meaning: Relic from an older resting state epoch strategy. 
    </details>
    </li>
    
  <li>
    <details>
      <summary>Head motion</summary>
      Where: fcp_1_TaskEpoching
      <br>
      Meaning: Specifics for initial handling of head motion
      <br><ul><li>
        <details>
        <summary>Threshold</summary>
        Where: fcp_1_TaskEpoching in head motion correction
        <br>
        Meaning: Threshold for which to reject trials with head motion 
        </details>
      </li></ul>
    </details>
    </li></ul>
</details>
