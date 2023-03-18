# Web users activities analytics

<br>Here is presented jupyter notebook with analysis of users activities. We want to indicate users by thier activities. 
This can be used for detecting a mailbox cracker by their behavior </br>

Project steps:
* goals and objectives of the project;
* description of the initial and processed data and features;
* primary feature analysis;
* primary visual feature analysis;
* description of data preprocessing;
* cross-validation, building validation and learning curves;
* evaluation of the model with a description of the selected metric;
* conclusions;

<br>The model gives good quality when distinguishing one user from another. It can be tried on, as mentioned above, to identify broken accounts. 
But for this I need to have a history of customer requests. There is also a flaw in the further training of the model. 
We made calculations on a limited number of users, and on a huge one it is still difficult for me to predict the result of accuracy. 
In addition, I am stubborn during the operation of the model itself, because a cracker can do his bad job, and the model will respond only after some time.</br>
