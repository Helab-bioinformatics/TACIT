import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import os
import matplotlib.pyplot as plt
import time

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

data_path = '/gpfs/share/home/t2406041101/xubin/20240831_new/train_set'
data_path2 = '/gpfs/share/home/t2406041101/xubin/20240831_new/target_set'
data_path3 = '/gpfs/share/home/t2406041101/xubin/20240831_new/train_set_1'
data_path4 = '/gpfs/share/home/t2406041101/xubin/20240831_new/target_set_1'
data_path5 = '/gpfs/share/home/t2406041101/xubin/20240831_new/train_set_2'
data_path6 = '/gpfs/share/home/t2406041101/xubin/20240831_new/target_set_2'

def load_data(input_paths, target_paths):
    input_data, target_data = [], []
    for input_file, target_file in zip(input_paths, target_paths):
        input_array = np.loadtxt(input_file)
        target_array = np.loadtxt(target_file)
        input_data.append(input_array)
        target_data.append(target_array)
    return torch.tensor(np.array(input_data), dtype=torch.float32).to(device), torch.tensor(np.array(target_data), dtype=torch.float32).to(device)

train_inputs = [os.path.join(data_path, f'input_matrix_{i}.txt') for i in range(1, 128)]
train_targets = [os.path.join(data_path2, f'output_matrix_{i}.txt') for i in range(1, 128)]
test_inputs1 = [os.path.join(data_path3, f'input_matrix_{i}.txt') for i in range(1, 32)]
test_targets1 = [os.path.join(data_path4, f'output_matrix_{i}.txt') for i in range(1, 32)]
test_inputs2 = [os.path.join(data_path5, f'input_matrix_{i}.txt') for i in range(1, 6)]
test_targets2 = [os.path.join(data_path6, f'output_matrix_{i}.txt') for i in range(1, 6)]

input_data, target_data = load_data(train_inputs, train_targets)
input_data1, target_data1 = load_data(test_inputs1, test_targets1)
input_data2, target_data2 = load_data(test_inputs2, test_targets2)

class CNNModel(nn.Module):
    def __init__(self):
        super(CNNModel, self).__init__()
        self.conv1 = nn.Conv1d(in_channels=3, out_channels=32, kernel_size=3, padding=1)
        self.pool = nn.MaxPool1d(kernel_size=2, stride=2)
        self.conv2 = nn.Conv1d(in_channels=32, out_channels=64, kernel_size=3, padding=1)
        self.fc1 = nn.Linear(64 * (273119 // 4), 512)  
        self.fc2 = nn.Linear(512, 256)
        self.fc3 = nn.Linear(256, 128)
        self.fc4 = nn.Linear(128, 64)
        self.fc5 = nn.Linear(64, 32)
        self.fc6 = nn.Linear(32, 3 * 273119)
        self.leaky_relu = nn.LeakyReLU(0.1)
        self.dropout = nn.Dropout(p=0.1)

    def forward(self, x):
        x = self.leaky_relu(self.conv1(x))
        x = self.pool(x)
        x = self.leaky_relu(self.conv2(x))
        x = self.pool(x)
        x = x.view(x.size(0), -1)  
        x = self.dropout(self.leaky_relu(self.fc1(x)))
        x = self.dropout(self.leaky_relu(self.fc2(x)))
        x = self.dropout(self.leaky_relu(self.fc3(x)))
        x = self.dropout(self.leaky_relu(self.fc4(x)))
        x = self.dropout(self.leaky_relu(self.fc5(x)))
        x = self.fc6(x)
        x = x.view(-1, 3, 273119)
        return x

model = CNNModel().to(device)

criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.005)


class EarlyStopping:
    def __init__(self, patience=200, delta=0.01):
        self.patience = patience
        self.delta = delta
        self.counter = 0
        self.best_loss = None
        self.early_stop = False

    def __call__(self, loss):
        if self.best_loss is None:
            self.best_loss = loss
        elif loss < self.best_loss - self.delta:
            self.best_loss = loss
            self.counter = 0
        else:
            self.counter += 1
            if self.counter >= self.patience:
                self.early_stop = True

early_stopping = EarlyStopping(patience=200)

losses = []
correlations = []
test_losses = []
test_correlations = []
test2_losses = []
test2_correlations = []


log_file = '20240418_training_log_cnn.txt'


total_start_time = time.time()


num_epochs = 1000
with open(log_file, 'w') as f:
    for epoch in range(num_epochs):
        start_time = time.time()  
        
        model.train()
        outputs = model(input_data)
        loss = criterion(outputs, target_data)
        
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
        end_time = time.time()  
        

        losses.append(loss.item())
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


        epoch_time = end_time - start_time
        gpu_memory_allocated = torch.cuda.memory_allocated(device) / (1024 ** 2) 
        gpu_max_memory_allocated = torch.cuda.max_memory_allocated(device) / (1024 ** 2)  

        f.write(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}, Correlation: {correlation:.4f}, '
                f'Test Loss: {test_loss.item():.4f}, Test Correlation: {test_correlation:.4f}, '
                f'Test embryos Loss: {test2_loss.item():.4f}, Test embryos Correlation: {test2_correlation:.4f}, '
                f'Time: {epoch_time:.2f}s, GPU Memory: {gpu_memory_allocated:.2f}MB, Max GPU Memory: {gpu_max_memory_allocated:.2f}MB\n')

        if (epoch+1) % 1 == 0:
            print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}, Correlation: {correlation:.4f}, '
                  f'Test Loss: {test_loss.item():.4f}, Test Correlation: {test_correlation:.4f}, '
                  f'Test embryos Loss: {test2_loss.item():.4f}, Test embryos Correlation: {test2_correlation:.4f}, '
                  f'Time: {epoch_time:.2f}s, GPU Memory: {gpu_memory_allocated:.2f}MB, Max GPU Memory: {gpu_max_memory_allocated:.2f}MB')

        # Early Stopping
        early_stopping(loss.item())
        if early_stopping.early_stop:
            print("Early stopping triggered")
            break


total_end_time = time.time()
total_training_time = total_end_time - total_start_time


with open('training_summary_cnn.txt', 'w') as summary_file:
    summary_file.write(f'Total Training Time: {total_training_time:.2f} seconds\n')
    summary_file.write(f'Max GPU Memory Allocated: {gpu_max_memory_allocated:.2f} MB\n')


torch.save(model, '20240809_save_cnn.pt')


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
plt.savefig('20240809_training_plots_cnn.svg')
plt.show()
