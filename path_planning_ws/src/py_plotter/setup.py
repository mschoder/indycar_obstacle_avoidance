from setuptools import setup

package_name = 'py_plotter'

setup(
    name=package_name,
    version='0.0.0',
    packages=[package_name],
    data_files=[
        ('share/ament_index/resource_index/packages',
            ['resource/' + package_name]),
        ('share/' + package_name, ['package.xml']),
    ],
    install_requires=['setuptools'],
    zip_safe=True,
    maintainer='mschoder',
    maintainer_email='schoder.m@gmail.com',
    description='plot path planning in rt',
    license='Apache License 2.0',
    tests_require=['pytest'],
    entry_points={
        'console_scripts': [
            'test_plotter = py_plotter.test_plotter:main'
        ],
    },
)
