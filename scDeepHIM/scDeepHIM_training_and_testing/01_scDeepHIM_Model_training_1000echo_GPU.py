import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import os
import matplotlib.pyplot as plt
import time  

# data loader
data_path = '/path1/input_training'
data_path2 = '/path2/output_training'
data_path3 = '/path3/input_testing'
data_path4 = '/path4/output_testing'
data_path5 = '/path5/input_validation'
data_path6 = '/path6/output_validation'

input_data = []
target_data = []
for i in range(1, 128):  
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

input_data1 = []
target_data1 = []
for i in range(1, 32):  
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

input_data2 = []
target_data2 = []
for i in range(1, 6):  
    input_file2 = os.path.join(data_path5, f'input_matrix_{i}.txt')
    target_file2 = os.path.join(data_path6, f'output_matrix_{i}.txt')
    
    input_array2 = np.loadtxt(input_file2)
    target_array2 = np.loadtxt(target_file2)
        
    input_data2.append(input_array2)
    target_data2.append(target_array2)

input_data2 = np.array(input_data2)
target_data2 = np.array(target_data2)

input_data2 = torch.tensor(input_data2, dtype=torch.float32)
target_data2 = torch.tensor(target_data2, dtype=torch.float32)


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


device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f'Using device: {device}')


model = LinearModel().to(device)


input_data = input_data.to(device)
target_data = target_data.to(device)
input_data1 = input_data1.to(device)
target_data1 = target_data1.to(device)
input_data2 = input_data2.to(device)
target_data2 = target_data2.to(device)


criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.005)

losses = []
learning_rates = []
correlations = []
test_losses = []
test_correlations = []
test2_losses = []
test2_correlations = []


start_time = time.time()


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
        output_array = outputs.detach().cpu().numpy().flatten()
        target_array = target_data.cpu().numpy().flatten()
        correlation = np.corrcoef(output_array, target_array)[0, 1]
        correlations.append(correlation)


    with torch.no_grad():
        model.eval()
        test_outputs = model(input_data1)
        test_loss = criterion(test_outputs, target_data1)
        test_output_array = test_outputs.cpu().numpy().flatten()
        test_target_array = target_data1.cpu().numpy().flatten()
        test_correlation = np.corrcoef(test_output_array, test_target_array)[0, 1]
        test_correlations.append(test_correlation)
    test_losses.append(test_loss.item())

    with torch.no_grad():
        model.eval()
        test_outputs2 = model(input_data2)
        test2_loss = criterion(test_outputs2, target_data2)
        test_output_array2 = test_outputs2.cpu().numpy().flatten()
        test_target_array2 = target_data2.cpu().numpy().flatten()
        test2_correlation = np.corrcoef(test_output_array2, test_target_array2)[0, 1]
        test2_correlations.append(test2_correlation)
    test2_losses.append(test2_loss.item())


    if torch.cuda.is_available():
        gpu_mem_allocated = torch.cuda.memory_allocated(device) / (1024 ** 2)  
        gpu_mem_reserved = torch.cuda.memory_reserved(device) / (1024 ** 2)  
        print(f'GPU Memory Allocated: {gpu_mem_allocated:.2f} MB, Reserved: {gpu_mem_reserved:.2f} MB')

    if (epoch+1) % 1 == 0:
        print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}, LearningRate: {current_lr:.6f}, Correlation: {correlation:.4f}, Test Loss: {test_loss.item():.4f}, Test Correlation: {test_correlation:.4f}, Test embryos Loss: {test2_loss.item():.4f}, Test embryos Correlation: {test2_correlation:.4f}')


end_time = time.time()
total_time = end_time - start_time
print(f'Total training time: {total_time:.2f} seconds')

torch.save(model,'Training_plots_DNN_1000echo.pt')
fig, axs = plt.subplots(2, 1, figsize=(8, 12))


axs[0].plot(losses, label='Training Loss')
axs[0].plot(test_losses, label='Testing Loss')
axs[0].plot(test2_losses, label='Testing embryos Loss')
axs[0].set_xlabel('Epoch')
axs[0].set_ylabel('Loss')
axs[0].set_title('Training Loss')
axs[0].legend()

axs[1].plot(correlations, label='Correlation')
axs[1].plot(test_correlations, label='Test_Correlation')
axs[1].plot(test2_correlations, label='Testing embryos Correlations')
axs[1].set_xlabel('Epoch')
axs[1].set_ylabel('Correlation')
axs[1].set_title('Correlation')
axs[1].legend()

plt.tight_layout()
plt.savefig('Training_plots_DNN_1000echo.svg')  
plt.show()
