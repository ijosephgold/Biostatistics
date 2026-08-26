#Intro to R

#This is an R script 
#Notice the hashtags before each line
#This tells the program to ignore these lines
#They are used for writing notes rather than code
#we'll get to code in a moment

#First, get familiar with the four windows:

#This is the script window
#Scripts are a useful way to write code and save it

#The window below this one is the console
#The console is actually the R program
#Once you write code, you'll notice that it runs in the console
#The output of many functions will be returned in the console

#The window at top right is the environment
#All data, values, and objects created in the code will be listed there
#It's a useful way to know what is going on with your code

#The window at bottom right shows the program files, and any graphics created
#It also shows any packages that have been installed

#Let's do a little basic coding
#The most basic application of R is to use it as a calculator
#Create some basic objects
#Notice we won't use hashtags
#create a variable x given a value of 1
x <- 1 #the "<-" syntax is used to create objects (it's kind of like "=")
x #call x by clicking Run in the top right
#Notice it pops up in the console and in the environment
#you can also just enter x into the console and hit enter
#create a variable y given a value of 2
y <- 2
#add x + y and call it z
z <- x + y
z #call z
#multiple x*y
xy <- x*y
xy
#divide xy by z
xy/z <- xy/z
#uh oh -- notice there is an error
#errors happen when there is improper syntax (among other issues)
#the slash in the object identifier is the problem
xyz <- xy/z
xyz #now it works
#Note that we could have done all of this in the console rather than in the script

#Let's build a matrix of random numbers
matrix <- matrix(c(1, 2, 3, 4, 5, 10, 20, 30, 40, 50), nrow=5, ncol=2)
matrix
#matrix() is a function
#We'll use lots of functions in the R language
#c() is a function that concatenates the specified values
#nrow and ncol specifies the number of rows and columns
#We could switch them around
matrix <- matrix(c(1, 2, 3, 4, 5, 10, 20, 30, 40, 50), nrow=2, ncol=5)
matrix
#Notice that we have redefined the object called matrix
#If we wanted to keep the original one, we could have named the new one something else
#If these numbers meant something, we could name them
row.names(matrix) <- c("Group 1", "Group 2") 
#Notice the names have to be put in "quotes" because they are characters
matrix
colnames(matrix) <- c("Alpha", "Beta", "Gamma", "Delta", "Epsilon")
matrix
#Creating a matrix is a useful tool if the data is not on a spreadsheet
#It's a way to create a data frame directly in R
#look at the structure of the matrix
str(matrix)
#convert it into a data frame
matrix <- as.data.frame(matrix)
str(matrix)
#Now we can call variables from it and make summary statistics
#For example
Alpha <- matrix$Alpha #the $ operator calls the Alpha column from the matrix data
#Calculate Alpha's mean
MeanAlpha <- mean(Alpha)
MeanAlpha
#Could also do that this way
MeanAlpha <- mean(matrix$Alpha)
#Calculate Alpha's standard deviation
StdevAlpha <- sd(matrix$Alpha)
StdevAlpha
#If we had data in a spreadsheet, we could load it into R
#But to do that, R needs to be able to find the file
#The best way to do this is to start an R project
#Click on the tab in the far upper right of R studio that says "Project (None)"
#Select "New Project" from the dropdown menu
#You'll need to name the project "Directory Name"
#Specify the directory location (pick a folder somewhere on your computer)
#This script will go away but you can reopen it in the new project
#When you get into the new project, go to File -> Recent Files -> Intro to R 
#Also open "R script template" in the new project
#Click create project
