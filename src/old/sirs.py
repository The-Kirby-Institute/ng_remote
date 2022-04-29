import random
import matplotlib.pyplot as plt


#%%###########################################################################
##  SETUP THE PERSON CLASS
##############################################################################


# MAKE THE PERSON CLASS
class person:
    id: int
    age: int
    sex: str
    location: int
    compartment: str
    compartment_t0: float
    compartment_dt: float

    # Define attributes of the person
    def __init__(self, id, age, sex, location, compartment, compartment_t0, compartment_dt = 7):
        self.id = id
        self.age = age
        self.sex = sex
        self.location = location
        self.compartment = compartment
        self.compartment_t0 = compartment_t0
        self.compartment_dt = compartment_dt

    # Set what to print when somebody asks about this object
    def __str__(self):
        return f"id: {self.id}\n age: {self.age}\n sex: {self.sex}\n location: {self.location}\n compartment:   {self.compartment}\n"

    # Instance method (a function which can be run)
    def description(self):
        return f"Person {self.id} is a {self.age} year old {self.sex} located in {self.location}"
    
    # Define the probability of the person causing a new infection
    def transmission_probability(self):
        return 0.1
    
    # Define the probability of the person causing a new infection
    def get_infectee(self, other_people):
        return random.choice(other_people)
    
    # Define method for calculating the duration of infection
    def duration_infectious(self):
        return 7
    
    # Define method for calculating the duration of immunity
    def duration_removed(self):
        return 30



#%%###########################################################################
##  INITILISE A LIST OF PEOPLE
##############################################################################
n_people = 1000
people = []
for i in range(0, n_people):
        
    # Determine the age of the individual
    age = random.uniform(0, 100)
    
    
    # Determine the sex of this individual
    if random.random() > 0.5:
        sex = "Male"
    else:
        sex = "Female"
        
        
    # Decide if they are infectious or not
    if random.random() > 0.01:
        compartment = "S"
    else:
        compartment = "I"
        
    
    # Make the individual
    p = person(id = i, 
               age = random.uniform(0, 100), 
               sex = sex, 
               location = 1, 
               compartment = compartment, 
               compartment_t0 = 0)
    p.compartment_dt = p.duration_infectious()
    
    
    # Append it to the list of people
    people.append(p)
    


#%%###########################################################################
##  SIMULATE AN INFECTION PROCESS
##############################################################################
t = range(0, 1000, 1)


# Initilise infection status for everyone
I = []
for p in people:
    if p.compartment == "I":
        I.append(p.id)


# Initilise S and R status for everyone else
S = [p.id for p in people if p.id not in I]
R = []
        

# Iterate over each time point
for t in t:
        
        
    # Work out if any infection events occur
    I_next = []
    for i in I:
        if len(S) > 0:
            if random.random() < people[i].transmission_probability():
                
                # Seed the infection
                I_new = p.get_infectee(S)
                people[I_new].compartment = "I"
                people[I_new].compartment_t0 = t
                people[I_new].compartment_dt = people[I_new].duration_infectious()
                
                # Update the compartment lists
                S.remove(I_new)
                I_next.append(I_new)
            
            
    # Update list of infectious individuals
    I = I + I_next
    
    
    # Determine if anybody is no longer infecious
    for i in I:
        if t > (people[i].compartment_t0 + people[i].compartment_dt):
            
            # Update the persons compartment status
            people[i].compartment = "R"
            people[I_new].compartment_t0 = t
            people[I_new].compartment_dt = people[I_new].duration_removed()
            
            # Update compartment lists
            I.remove(i)
            R.append(i)
            
            
    # Determine if anybody is no longer removed
    for r in R:
        if t > (people[r].compartment_t0 + people[r].compartment_dt):
            
            # Update the persons compartment status
            people[r].compartment = "S"
            people[r].compartment_t0 = t
            people[r].compartment_dt = 0
            
            # Update compartment lists
            R.remove(r)
            S.append(r)
            
            
    plt.scatter([t, t, t], [len(S), len(I), len(R)], c=["r", "g", "b"])
    
    
    
                
            






#%%
















