     
import pandas as pd


df = pd.DataFrame(
    data_dict
    )
     
#read data
path = r"C:\Users\MRahim\Downloads\2026Spring_BE_7910_Course_Materials\Assignment\Data\titanic"
train_df = pd.read_csv("train.csv")


#subset

#only columns
age = train_df["Age"]

custom_df = train_df[["Age","Fare","Pclass"]]


#iloc
train_df.iloc[0:11,2:4]


#loc
train_df.loc[0:10,["Name","Pclass","Age"]]


train_df.set_index("PassengerId")

train_df.set_index("PassengerId",inplace=True)
train_df.reset_index(inplace=True)


#check condition

train_df["Age"] >50
train_df["Fare"] > 80

df_age_50 = train_df[train_df["Age"] >50]


a = train_df["Age"] >50
b = train_df["Fare"] > 80

df_age_fare = train_df[(a) & (b)]


#using list comprehension

train_df["sample"] = 0


train_df["AgeGroup"] = ["Kid" if age<18 else "Senior" if age>60 else "Adult" for age in train_df["Age"]]

#using numpy
train_df["AgeGroup"] = np.where(
    train_df["Age"] < 18, "Child",
    np.where(train_df["Age"] > 60, "Senior", "Adult")
)


#column names
train_df.columns

#change column name
columns_name = train_df.columns.tolist()

columns_name[-2] = "Port Embarked"

train_df.columns = columns_name

#add new column

train_df["None"] = 0



#set index
train_df.set_index("Pclass", inplace=True)


#loops
for i in range(0,len(train_df)):
    train_df.iloc[i]
    break


for i in train_df.iterrows():
    break

#check na

train_df.isna().sum()

train_df["Age"].isna().sum()


#drop na

sample_df = train_df.dropna(axis=1, how = "any")


#fill na
df1 = train_df.fillna(0)

mean_age = 35

train_df["Age"].fillna(35)


#merge dataframes

df1 = train_df[["PassengerId", "Age"]]

df2 = train_df[["PassengerId", "Fare"]]

df3 = df1.merge(df2,left_on="PassengerId", right_on ="PassengerId" )


#concat

df4 = pd.concat([df1,df2],axis=1)



#write the csv

df3.to_csv(path+"/train_process.csv")


