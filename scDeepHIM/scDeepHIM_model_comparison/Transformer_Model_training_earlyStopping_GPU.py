import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import os
import matplotlib.pyplot as plt
import time
import torch.cuda.amp as amp  


device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


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


input_data = torch.tensor(input_data, dtype=torch.float32).to(device)
target_data = torch.tensor(target_data, dtype=torch.float32).to(device)


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

input_data1 = torch.tensor(input_data1, dtype=torch.float32).to(device)
target_data1 = torch.tensor(target_data1, dtype=torch.float32).to(device)


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

input_data2 = torch.tensor(input_data2, dtype=torch.float32).to(device)
target_data2 = torch.tensor(target_data2, dtype=torch.float32).to(device)


losses = []
test_losses = []
test2_losses = []
correlations = []
test_correlations = []
test2_correlations = []

class TransformerModel(nn.Module):
    def __init__(self, input_dim=3, seq_length=273119, num_heads=2, num_layers=1, dim_feedforward=32, dropout=0.1):
        super(TransformerModel, self).__init__()
        self.input_linear = nn.Linear(input_dim, dim_feedforward)
        self.positional_encoding = nn.Parameter(torch.zeros(1, seq_length, dim_feedforward))
        encoder_layer = nn.TransformerEncoderLayer(d_model=dim_feedforward, nhead=num_heads, 
                                                   dim_feedforward=dim_feedforward, dropout=dropout, batch_first=True)
        self.transformer_encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        self.output_linear = nn.Linear(dim_feedforward, input_dim)
        self.leaky_relu = nn.LeakyReLU(0.1)
        self.dropout = nn.Dropout(p=dropout)

    def forward(self, x):
        batch_size, channels, seq_length = x.size()  
        x = x.permute(0, 2, 1)  
        x = self.input_linear(x)
        

        pos_encoding = self.positional_encoding[:, :seq_length, :]
        
        x += pos_encoding
        x = self.dropout(x)
        x = self.transformer_encoder(x)
        x = self.output_linear(x)
        x = x.permute(0, 2, 1)  
        x = self.leaky_relu(x)
        return x

model = TransformerModel().to(device)


criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.005)


scaler = amp.GradScaler()


log_file = 'Training_log.txt'


total_start_time = time.time()


patience = 50
delta = 0.01
best_loss = float('inf')
counter = 0
early_stop = False


num_epochs = 1000
chunk_size = 10000  

with open(log_file, 'w') as f:
    for epoch in range(num_epochs):
        if early_stop:
            print(f'Early stopping at epoch {epoch+1}')
            break
        
        start_time = time.time()  
        
        model.train()
        optimizer.zero_grad()
        
        total_loss = 0
        for start in range(0, input_data.size(2), chunk_size):
            end = min(start + chunk_size, input_data.size(2))
            input_chunk = input_data[:, :, start:end]
            target_chunk = target_data[:, :, start:end]
            

            with amp.autocast():
                outputs = model(input_chunk)
                loss = criterion(outputs, target_chunk)
            
            scaler.scale(loss).backward()
            total_loss += loss.item()
        
        scaler.step(optimizer)
        scaler.update()
        
        avg_loss = total_loss / (input_data.size(2) // chunk_size)  
        losses.append(avg_loss)


        with torch.no_grad():
            train_outputs = model(input_data)
            train_output_array = train_outputs.cpu().numpy().flatten()
            train_target_array = target_data.cpu().numpy().flatten()
            train_correlation = np.corrcoef(train_output_array, train_target_array)[0, 1]
            correlations.append(train_correlation)


        with torch.no_grad():
            test_loss1 = 0
            for start in range(0, input_data1.size(2), chunk_size):
                end = min(start + chunk_size, input_data1.size(2))
                input_chunk = input_data1[:, :, start:end]
                target_chunk = target_data1[:, :, start:end]
                
                test_outputs = model(input_chunk)
                loss = criterion(test_outputs, target_chunk)
                test_loss1 += loss.item()
            
            test_loss1 /= (input_data1.size(2) // chunk_size)
            test_losses.append(test_loss1)
            

            output_array = test_outputs.cpu().numpy().flatten()
            target_array = target_chunk.cpu().numpy().flatten()
            correlation1 = np.corrcoef(output_array, target_array)[0, 1]
            test_correlations.append(correlation1)
        

        with torch.no_grad():
            test_loss2 = 0
            for start in range(0, input_data2.size(2), chunk_size):
                end = min(start + chunk_size, input_data2.size(2))
                input_chunk = input_data2[:, :, start:end]
                target_chunk = target_data2[:, :, start:end]
                
                test_outputs = model(input_chunk)
                loss = criterion(test_outputs, target_chunk)
                test_loss2 += loss.item()
            
            test_loss2 /= (input_data2.size(2) // chunk_size)
            test2_losses.append(test_loss2)
            

            output_array2 = test_outputs.cpu().numpy().flatten()
            target_array2 = target_chunk.cpu().numpy().flatten()
            correlation2 = np.corrcoef(output_array2, target_array2)[0, 1]
            test2_correlations.append(correlation2)


        end_time = time.time()  
        epoch_time = end_time - start_time
        gpu_memory_allocated = torch.cuda.memory_allocated(device) / (1024 ** 2)  
        gpu_max_memory_allocated = torch.cuda.max_memory_allocated(device) / (1024 ** 2)  

        f.write(f'Epoch [{epoch+1}/{num_epochs}], Loss: {avg_loss:.4f}, Train Correlation: {train_correlation:.4f}, '
                f'Test Loss1: {test_loss1:.4f}, Test Correlation1: {correlation1:.4f}, '
                f'Test Loss2: {test_loss2:.4f}, Test Correlation2: {correlation2:.4f}, '
                f'Time: {epoch_time:.2f}s, GPU Memory: {gpu_memory_allocated:.2f}MB, Max GPU Memory: {gpu_max_memory_allocated:.2f}MB\n')

        if (epoch+1) % 1 == 0:
            print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {avg_loss:.4f}, Train Correlation: {train_correlation:.4f}, '
                  f'Test Loss1: {test_loss1:.4f}, Test Correlation1: {correlation1:.4f}, '
                  f'Test Loss2: {test_loss2:.4f}, Test Correlation2: {correlation2:.4f}, '
                  f'Time: {epoch_time:.2f}s, GPU Memory: {gpu_memory_allocated:.2f}MB, Max GPU Memory: {gpu_max_memory_allocated:.2f}MB')


        if avg_loss < best_loss - delta:
            best_loss = avg_loss
            counter = 0
        else:
            counter += 1
            if counter >= patience:
                early_stop = True


total_end_time = time.time()
total_training_time = total_end_time - total_start_time


with open('training_summary_transformer.txt', 'w') as summary_file:
    summary_file.write(f'Total Training Time: {total_training_time:.2f} seconds\n')
    summary_file.write(f'Max GPU Memory Allocated: {gpu_max_memory_allocated:.2f} MB\n')


torch.save(model, 'Save_transformer.pt')


fig, axs = plt.subplots(2, 1, figsize=(8, 12))


axs[0].plot(losses, label='Training Loss')
axs[0].plot(test_losses, label='Testing Loss1')
axs[0].plot(test2_losses, label='Testing Loss2')
axs[0].set_xlabel('Epoch')
axs[0].set_ylabel('Loss')
axs[0].set_title('Training and Testing Loss')
axs[0].legend()


axs[1].plot(correlations, label='Training Correlation')
axs[1].plot(test_correlations, label='Test Correlation1')
axs[1].plot(test2_correlations, label='Test Correlation2')
axs[1].set_xlabel('Epoch')
axs[1].set_ylabel('Correlation')
axs[1].set_title('Correlation')
axs[1].legend()

plt.tight_layout()
plt.savefig('Training_plots_transformer.svg')
plt.show()
