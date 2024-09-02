import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import os
import matplotlib.pyplot as plt
import time  # 用于计算时间

# 设置数据目录
data_path = '/gpfs/share/home/t2406041101/xubin/20240831_new/train_set'
data_path2 = '/gpfs/share/home/t2406041101/xubin/20240831_new/target_set'
data_path3 = '/gpfs/share/home/t2406041101/xubin/20240831_new/train_set_1'
data_path4 = '/gpfs/share/home/t2406041101/xubin/20240831_new/target_set_1'
data_path5 = '/gpfs/share/home/t2406041101/xubin/20240831_new/train_set_2'
data_path6 = '/gpfs/share/home/t2406041101/xubin/20240831_new/target_set_2'

# 加载数据(traing group)
input_data = []
target_data = []
for i in range(1, 128):  # 加载3套数据
    input_file = os.path.join(data_path, f'input_matrix_{i}.txt')
    target_file = os.path.join(data_path2, f'output_matrix_{i}.txt')
    
    # 读取txt文件并将数据转换为numpy数组
    input_array = np.loadtxt(input_file)
    target_array = np.loadtxt(target_file)
        
    input_data.append(input_array)
    target_data.append(target_array)

input_data = np.array(input_data)
target_data = np.array(target_data)

# 将数据转换为张量
input_data = torch.tensor(input_data, dtype=torch.float32)
target_data = torch.tensor(target_data, dtype=torch.float32)

# 加载数据(testing group)
input_data1 = []
target_data1 = []
for i in range(1, 32):  # 加载3套数据
    input_file1 = os.path.join(data_path3, f'input_matrix_{i}.txt')
    target_file1 = os.path.join(data_path4, f'output_matrix_{i}.txt')
    
    # 读取txt文件并将数据转换为numpy数组
    input_array1 = np.loadtxt(input_file1)
    target_array1 = np.loadtxt(target_file1)
        
    input_data1.append(input_array1)
    target_data1.append(target_array1)

input_data1 = np.array(input_data1)
target_data1 = np.array(target_data1)

# 将数据转换为张量
input_data1 = torch.tensor(input_data1, dtype=torch.float32)
target_data1 = torch.tensor(target_data1, dtype=torch.float32)

# 加载数据(testing group2)
input_data2 = []
target_data2 = []
for i in range(1, 6):  # 加载3套数据
    input_file2 = os.path.join(data_path5, f'input_matrix_{i}.txt')
    target_file2 = os.path.join(data_path6, f'output_matrix_{i}.txt')
    
    # 读取txt文件并将数据转换为numpy数组
    input_array2 = np.loadtxt(input_file2)
    target_array2 = np.loadtxt(target_file2)
        
    input_data2.append(input_array2)
    target_data2.append(target_array2)

input_data2 = np.array(input_data2)
target_data2 = np.array(target_data2)

# 将数据转换为张量
input_data2 = torch.tensor(input_data2, dtype=torch.float32)
target_data2 = torch.tensor(target_data2, dtype=torch.float32)


# 定义卷积神经网络模型
class LinearModel(nn.Module):
    def __init__(self):
        super(LinearModel, self).__init__()
        self.fc1 = nn.Linear(3 * 273119, 512)  # 增加隐藏层的神经元数量
        self.bn1 = nn.BatchNorm1d(512)  # 添加批归一化层
        self.dropout1 = nn.Dropout(p=0.1)  # 
        self.fc2 = nn.Linear(512, 256)  # 增加隐藏层的神经元数量
        self.bn2 = nn.BatchNorm1d(256)  # 添加批归一化层
        self.dropout2 = nn.Dropout(p=0.1)  # 
        self.fc3 = nn.Linear(256, 128)  # 增加隐藏层的神经元数量
        self.bn3 = nn.BatchNorm1d(128)
        self.dropout3 = nn.Dropout(p=0.1)  # 
        self.fc4 = nn.Linear(128, 64)  # 输出层
        self.bn4 = nn.BatchNorm1d(64)
        self.dropout4 = nn.Dropout(p=0.1)  # 
        self.fc5 = nn.Linear(64, 32)  # 输出层
        self.bn5 = nn.BatchNorm1d(32)
        self.dropout5 = nn.Dropout(p=0.1)  # 
        self.fc6 = nn.Linear(32, 3 * 273119)  # 输出层
        self.bn6 = nn.BatchNorm1d(3 * 273119)
        self.leaky_relu = nn.LeakyReLU(0.1)  # 使用Leaky ReLU激活函数

    def forward(self, x):
        x = x.view(-1, 3 * 273119)  # 展平输入数据
        x = self.leaky_relu(self.bn1(self.fc1(x)))  # 使用Leaky ReLU激活函数
        x = self.dropout1(x)  # 在全连接层后添加dropout
        x = self.leaky_relu(self.bn2(self.fc2(x)))  # 使用Leaky ReLU激活函数
        x = self.dropout2(x)  # 在全连接层后添加dropout
        x = self.leaky_relu(self.bn3(self.fc3(x)))  # 使用Leaky ReLU激活函数
        x = self.dropout3(x)  # 在全连接层后添加dropout
        x = self.leaky_relu(self.bn4(self.fc4(x)))  # 使用Leaky ReLU激活函数
        x = self.dropout4(x)  # 在全连接层后添加dropout
        x = self.leaky_relu(self.bn5(self.fc5(x)))  # 使用Leaky ReLU激活函数
        x = self.dropout5(x)  # 在全连接层后添加dropout
        x = self.leaky_relu(self.bn6(self.fc6(x)))  # 使用Leaky ReLU激活函数
        x = x.view(-1, 3, 273119)  # 调整输出
        return x

# 设置设备为 GPU（如果可用）
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f'Using device: {device}')

# 将模型移到 GPU
model = LinearModel().to(device)

# 将数据移到 GPU
input_data = input_data.to(device)
target_data = target_data.to(device)
input_data1 = input_data1.to(device)
target_data1 = target_data1.to(device)
input_data2 = input_data2.to(device)
target_data2 = target_data2.to(device)

# 定义损失函数和优化器
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.005)

losses = []
learning_rates = []
correlations = []
test_losses = []
test_correlations = []
test2_losses = []
test2_correlations = []

# 记录开始时间
start_time = time.time()

# 训练模型
num_epochs = 1000
for epoch in range(num_epochs):
    # 前向传播
    model.train()
    outputs = model(input_data)
    loss = criterion(outputs, target_data)
    
    # 反向传播和优化
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()

    losses.append(loss.item())

    current_lr = optimizer.param_groups[0]['lr']

    # 计算相关性时需要断开梯度追踪
    with torch.no_grad():
        output_array = outputs.detach().cpu().numpy().flatten()
        target_array = target_data.cpu().numpy().flatten()
        correlation = np.corrcoef(output_array, target_array)[0, 1]
        correlations.append(correlation)

    # 计算测试集的各参数
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

    # 打印 GPU 内存占用情况
    if torch.cuda.is_available():
        gpu_mem_allocated = torch.cuda.memory_allocated(device) / (1024 ** 2)  # 转换为 MB
        gpu_mem_reserved = torch.cuda.memory_reserved(device) / (1024 ** 2)  # 转换为 MB
        print(f'GPU Memory Allocated: {gpu_mem_allocated:.2f} MB, Reserved: {gpu_mem_reserved:.2f} MB')

    if (epoch+1) % 1 == 0:
        print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}, LearningRate: {current_lr:.6f}, Correlation: {correlation:.4f}, Test Loss: {test_loss.item():.4f}, Test Correlation: {test_correlation:.4f}, Test embryos Loss: {test2_loss.item():.4f}, Test embryos Correlation: {test2_correlation:.4f}')

# 记录结束时间
end_time = time.time()
total_time = end_time - start_time
print(f'Total training time: {total_time:.2f} seconds')

torch.save(model,'20240418_save_3layers_1000echo_0901.pt')
fig, axs = plt.subplots(2, 1, figsize=(8, 12))

# 可视化损失变化
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
plt.savefig('20240901_training_plots_3layers_1000echo.svg')  # 保存图形
plt.show()
