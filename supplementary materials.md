Dear reviewer,  
&ensp;&ensp;We are very sorry that the unclear expression in the first edition of the manuscript wasted your time and didn't clearly show you our method. Thank you very much for your question. These problems are related to the accuracy of our proposed method. Let's solemnly and carefully check our process, and carry out various relevant experiments to verify it, so as to show our method to readers better.  
&ensp;&ensp;For the first problem, we also found that Zhang Song's paper also uses the fitting method to obtain the ideal phase ("flexible and high accuracy method for uni directional structured light system calibration"). His fitting way is different from ours as follows:

&ensp;&ensp;&ensp;&ensp;1) The primary sources of phase error are different  
&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;His error mainly comes from the system's noise and the calibration plate's colour change, and we are gamma effect.  <br>   &ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;However, their error margin is not lower than ours.
![image](https://user-images.githubusercontent.com/61692154/168427980-fd16c3e5-ffad-403e-89a3-22a8c1c956e4.png)

&ensp;&ensp;&ensp;&ensp;2) The use of the ideal phase is different  
&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;They are mainly used for system calibration, and then use a new algorithm for reconstruction, which requires higher accuracy than our <br> &ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;method. During reconstruction, the phase is directly brought into the following equation to obtain the actual spatial coordinates of <br> &ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;the point.

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;![image](https://user-images.githubusercontent.com/61692154/168428174-c0eda65a-e01e-4046-a4b6-48708e05d398.png)

&ensp;&ensp;&ensp;&ensp;3) Different fitting methods and objects  
&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;Because they take it for calibration, they have to fit the whole plate (calibration plate) with different attitudes. We only need a small area <br>&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;in the centre of the plate, which is perpendicular to the optical axis of the projector. Zhang's fitting way follows as follows:

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;![image](https://user-images.githubusercontent.com/61692154/168426951-6ecd6cc5-b2c8-4b90-9cd0-4f76cb47198a.png)

&ensp;&ensp;Their final results also prove the feasibility of this method.

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;![image](https://user-images.githubusercontent.com/61692154/168428658-3dc74a64-511a-4ab8-910a-f662841d5326.png)

&ensp;&ensp;Finally, we apologize to you again for our first edition, and please don't mind us showing you the content of other articles in the form of screenshots. We wish you good health and success in your work.

<p align="right">WangJie</p>
<p align="right">2022.05.14</p>
