# JSON Config Paramters

This file is the documentation for the paramters found in the JSON config file. Click the arrow on each parameter to reveal a dropdown explanation of where the paramter is used in the code and what it means. 

<details>
<summary>Contact</summary>
<br>
Where: sendEmail function is called at the end of each fcp_# step, and contact is passed as a parameter.
<br>
Meaning: Email address to which to send pipeline’s progress updates (contained in square brackets, [ ])
</details>

<details>
<summary>Epoching</summary>
<br>
Where: fcp_1_TaskEpoching.
<br>
Meaning: Epoch the data into trials
  
  <details>
  <summary>Period</summary>
  <br>
  Where: fcp_1_RestingStateEpoching, line 97
  <br>
  Meaning: Indicates epoch length for epoching resting state data
  </details>
  
  <details>
  <summary>Total time</summary>
  <br>
  Where: Nowhere
  <br>
  Meaning: Relic from an older resting state epoch strategy. 
  </details>
  
  <details>
  <summary>Head motion</summary>
  <br>
  Where: fcp_1_TaskEpoching 
  <br>
  Meaning:Specifics for initial handling of head motion
    <details> 
    <summary>Threshold</summary>
    <br>
    Where: fcp_1_TaskEpoching in head motion correction
    <br>
    Meaning: Threshold for which to reject trials with head motion 
    </details>
    
  </details>
  
</details>
