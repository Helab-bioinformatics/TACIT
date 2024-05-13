import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import os
import matplotlib.pyplot as plt

# Setting the data path
data_path = '/data_path1'
data_path2 = '/data_path2'
data_path3 = '/data_path3'
data_path4 = '/data_path4'

# Loading data for traing group (3*273119 matrics)
input_data = []
target_data = []
for i in range(1, n):   # n denotes the number of datasets for training
    input_file = os.path.join(data_path, f'input_matrix_{i}.txt')
    target_file = os.path.join(data_path2, f'output_matrix_{i}.txt')
    
    input_array = np.loadtxt(input_file)
    target_array = np.loadtxt(target_file)
        
    input_data.append(input_array)
    target_data.append(target_array)

input_data = np.array(input_data)
target_data = np.array(target_data)

input_data = torch.tensor(input_data, dtype=torch.float32)
target_data = torch.tensor(target_data, dtype=torch.float32)

# Loading data for testing group (3*273119 matrics)
input_data1 = []
target_data1 = []
for i in range(1, i):  # i denotes the number of datasets for testing
    input_file1 = os.path.join(data_path3, f'input_matrix_{i}.txt')
    target_file1 = os.path.join(data_path4, f'output_matrix_{i}.txt')
    
    input_array1 = np.loadtxt(input_file1)
    target_array1 = np.loadtxt(target_file1)
        
    input_data1.append(input_array1)
    target_data1.append(target_array1)

input_data1 = np.array(input_data1)
target_data1 = np.array(target_data1)

input_data1 = torch.tensor(input_data1, dtype=torch.float32)
target_data1 = torch.tensor(target_data1, dtype=torch.float32)

# Defining the DNN model
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
        
model = LinearModel()
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.005)

losses = []
learning_rates = []
correlations = []
test_losses = []
test_correlations = []

num_epochs = 1000
for epoch in range(num_epochs):
    model.train()
    outputs = model(input_data)
    loss = criterion(outputs, target_data)
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()
    losses.append(loss.item())
    current_lr = optimizer.param_groups[0]['lr']

    with torch.no_grad():
        output_array = outputs.detach().numpy().flatten()
        target_array = target_data.numpy().flatten()
        correlation = np.corrcoef(output_array, target_array)[0, 1]
        correlations.append(correlation)

    with torch.no_grad():
        model.eval()
        test_outputs = model(input_data1)
        test_loss = criterion(test_outputs, target_data1)
        test_output_array = test_outputs.numpy().flatten()
        test_target_array = target_data1.numpy().flatten()
        test_correlation = np.corrcoef(test_output_array, test_target_array)[0, 1]
        test_correlations.append(test_correlation)
    test_losses.append(test_loss.item())

    if (epoch+1) % 1 == 0:
        print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}, LearningRate: {current_lr:.6f}, Correlation: {correlation:.4f}，Test Loss: {test_loss.item():.4f}, Test Correlation: {test_correlation:.4f}')


torch.save(model,'scDeepHIM_1000echo.pt')
fig, axs = plt.subplots(2, 1, figsize=(8, 12))

#Visualization of the training process
axs[0].plot(losses, label='Training Loss')
axs[0].plot(test_losses, label='Testing Loss')
axs[0].set_xlabel('Epoch')
axs[0].set_ylabel('Loss')
axs[0].set_title('Training Loss')
axs[0].legend()

axs[1].plot(correlations, label='Correlation')
axs[1].plot(test_correlations, label='Test_Correlation')
axs[1].set_xlabel('Epoch')
axs[1].set_ylabel('Correlation')
axs[1].set_title('Correlation')
axs[1].legend()

plt.tight_layout()
plt.savefig('scDeepHIM_1000echo.svg')  # 保存图形
plt.show()



