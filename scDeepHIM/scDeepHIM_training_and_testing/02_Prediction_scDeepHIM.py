import torch
import torch.nn as nn
import numpy as np
import os

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

class LinearModel(nn.Module):
    def __init__(self):
        super(LinearModel, self).__init__()
        self.fc1 = nn.Linear(3 * 273119, 512)  
        self.bn1 = nn.BatchNorm1d(512)  
        self.dropout1 = nn.Dropout(p=0.1)   
        self.fc2 = nn.Linear(512, 256)  
        self.bn2 = nn.BatchNorm1d(256)  
        self.dropout2 = nn.Dropout(p=0.1)   
        self.fc3 = nn.Linear(256, 128) 
        self.bn3 = nn.BatchNorm1d(128)
        self.dropout3 = nn.Dropout(p=0.1)   
        self.fc4 = nn.Linear(128, 64)  
        self.bn4 = nn.BatchNorm1d(64)
        self.dropout4 = nn.Dropout(p=0.1)  
        self.fc5 = nn.Linear(64, 32)  
        self.bn5 = nn.BatchNorm1d(32)
        self.dropout5 = nn.Dropout(p=0.1) 
        self.fc6 = nn.Linear(32, 3 * 273119)  
        self.bn6 = nn.BatchNorm1d(3 * 273119)
        self.leaky_relu = nn.LeakyReLU(0.1)  

    def forward(self, x):
        x = x.view(-1, 3 * 273119)  
        x = self.leaky_relu(self.bn1(self.fc1(x)))  
        x = self.dropout1(x)  
        x = self.leaky_relu(self.bn2(self.fc2(x)))  
        x = self.dropout2(x)  
        x = self.leaky_relu(self.bn3(self.fc3(x)))  
        x = self.dropout3(x)  
        x = self.leaky_relu(self.bn4(self.fc4(x)))  
        x = self.dropout4(x)  
        x = self.leaky_relu(self.bn5(self.fc5(x)))  
        x = self.dropout5(x)  
        x = self.leaky_relu(self.bn6(self.fc6(x)))  
        x = x.view(-1, 3, 273119)  
        return x

# Load model
model = torch.load('scDeepHIM_model_1000echo_GPU.pt')
model.eval()

# Load data
data_path = 'path1' #path for the data to be predicted (with H3K27ac, H3K27me3 and H3K9me3) 
input_data = []
for i in range(1, n):  #n represents the number of samples to be predicted
    input_file = os.path.join(data_path, f'input_matrix_{i}.txt')
    
    input_array = np.loadtxt(input_file)
        
    input_data.append(input_array)

input_data = np.array(input_data)

# Convert to PyTorch tensors and move to the same device as the model
input_data = torch.tensor(input_data, dtype=torch.float32).to(device)

# Make predictions
with torch.no_grad():
    output = model(input_data)

predictions = output.cpu().numpy()  # Move predictions to CPU if necessary
num_samples = predictions.shape[0]
for i in range(num_samples):
    prediction_i = predictions[i]
    np.savetxt(f'prediction_sample_{i+1}.txt', prediction_i)
