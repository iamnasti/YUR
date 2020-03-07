from peewee import *
import pandas
import sqlite3

DATABASE = 'base.db'

database = sqlite3.connect(DATABASE)


class BaseModel(Model):
    class Meta:
        database = database


class Training(BaseModel):
    Training_Id = IntegerField(primary_key=True)
    Target = CharField()
    Date = DateTimeField()


class Competition(BaseModel):
    Competition_Id = IntegerField(primary_key=True)
    Name = CharField()
    Date = DateTimeField()
    Duration = IntegerField()


class Coach(BaseModel):
    Coach_Id = IntegerField(primary_key=True)
    Name = CharField()
    Experience = IntegerField()


class Prize(BaseModel):
    Prize_Id = IntegerField(primary_key=True)
    Name = CharField()
    Title = CharField()


class Sportsman(BaseModel):
    Sportsman_Id = IntegerField(primary_key=True)
    Name = CharField()
    Sport_Name = CharField()
    Coach_Id = ForeignKeyField(Coach)


class Fan(BaseModel):
    Fan_Id = IntegerField(primary_key=True)
    Name = CharField()
    Club = CharField()
    Sportsman_Id = ForeignKeyField(Sportsman)


class M_Competition(BaseModel):
    id_sportsman = ForeignKeyField(Sportsman)
    id_competition = ForeignKeyField(Competition)


class M_Traning(BaseModel):
    id_sportsman = ForeignKeyField(Sportsman)
    id_traning = ForeignKeyField(Training)


csv_path_fan = 'data_fan.csv'
csv_path_coach = 'data_coach.csv'
csv_path_comp = 'data_comp.csv'
csv_path_m_competition = 'data_m_competition.csv'
csv_path_m_training = 'data_m_training.csv'
csv_path_prize = 'data_prize.csv'
csv_path_sportsman = 'data_sportsman.csv'
csv_path_traning = 'data_traning.csv'
pandas.read_csv(csv_path_traning).to_sql("Training", database, if_exists='append', index=False)
pandas.read_csv(csv_path_comp).to_sql("Competition", database, if_exists='append', index=False)
pandas.read_csv(csv_path_coach).to_sql("Coach", database, if_exists='append', index=False)
pandas.read_csv(csv_path_prize).to_sql("Prize", database, if_exists='append', index=False)
pandas.read_csv(csv_path_sportsman).to_sql("Sportsman", database, if_exists='append', index=False)
pandas.read_csv(csv_path_fan).to_sql("Fan", database, if_exists='append', index=False)
pandas.read_csv(csv_path_m_competition).to_sql("M_Competition", database, if_exists='append', index=False)
pandas.read_csv(csv_path_m_training).to_sql("M_Traning", database, if_exists='append', index=False)